import datajoint as dj
import numpy as np 
import random
import opexebo
from tqdm.auto import tqdm


# Load base schema
schema = dj.schema(dj.config['dj_imaging.database'])
schema.spawn_missing_classes()

# Load personal schema 
borderscore_schema = dj.schema('user_horsto_borderscore')
borderscore_schema.spawn_missing_classes()

from bvs.detect_fields import detect_fields
from bvs.bv_score import calc_bv_score

from helpers.utils import find_nearest, calc_signal_dict, calc_ratemap


@borderscore_schema
class ShuffledBVS(dj.Computed):
    definition = """
    # Shuffling table for boundary vector score (BVS)
    -> FilteredSpikes.proj(signal_dataset = 'dataset_name')
    -> Tracking.OpenField.proj(tracking_dataset = 'dataset_name')
    -> ShuffleParams
    -> BVField
    -> BVScoreFieldMethod
    ---
    number_shuffles            :    int                   # Total number of shuffles (can vary from expected number)
    shuffling_offsets          :    blob@imgstore         # Shuffling offsets
    """

    class BVS(dj.Part):
        definition = """
        # Shuffled bvs  
        -> master
        --- 
        bvs_99                 :  double         #  Information rate 99th percentile
        bvs_95                 :  double         #  Information rate 95th percentile
        bvs_shuffles           :  blob@imgstore  #  Individual shuffles information rate   
        """


    def make(self,key):

        ''' 
        This routine follows the basic layout of "SignalTracking" (see ratemaps.py in main imaging schema).
        It extracts sync and signal data and creates a vector of shuffling integers 
        that are used to 'roll' the frames sync signal around. 
        This maintains the underlying signal statistics, and destroys the correspondence
        between signal and tracking data.
        '''    

        bvs_field_method = (BVScoreFieldMethod & key).fetch1('bv_field_dect_method')
        if bvs_field_method not in ['opexebo','bvs']:
            raise NotImplementedError(f'Method "{bvs_field_method}" does not have a matching score calculation routine')


        shuffle_params = (ShuffleParams & key).fetch1()
        st_params = {}
        st_params['speed_cutoff_low'], st_params['speed_cutoff_high'], st_time_offset = (SignalTrackingParams & key).fetch1(
                                                                    'speed_cutoff_low', 'speed_cutoff_high', 'time_offset')

        spikes   = (FilteredSpikes.proj(signal_dataset='dataset_name', spikes='filtered_spikes')
                    & key).fetch1('spikes')
        tracking = (Tracking.OpenField * Tracking.proj(tracking_dataset='dataset_name')
                    & key).fetch1()

        center_y, center_plane = (Cell.Rois.proj(..., signal_dataset='dataset_name') & key).fetch1(
                    'center_y', 'center_plane')
        num_planes = (Tif.SI & (Session & key)).fetch1('num_scanning_depths')

        # Ratemap 
        occupancy_entry  = (Occupancy & key).fetch1()
        ratemap_params   = (MapParams & key).fetch1()
        field_params     = (FieldParams & key).fetch1()

        # BV Field method 
        bvs_field_params =  (BVFieldParams & key).fetch1()

        occupancy = np.ma.array(occupancy_entry['occupancy'], mask=occupancy_entry['mask_occ'])
        x_edges = occupancy_entry['x_edges'].copy()
        y_edges = occupancy_entry['y_edges'].copy()

        # Experiment Type
        experiment_type   = (Session & key).fetch1('experiment_type')


        # Retrieve Sync
        if (TrackingRaw & key).fetch1('sync'):
            sync_data_frames, sample_rate_frames = (MetaSession.Setup * Setup.Sync \
                                                  * Sync & 'generic_name = "frames_imaging"' & key).fetch1(\
                                                  'sync_data', 'sample_rate')
            sync_data_track = (MetaSession.Setup * Setup.Sync \
                                                 * Sync & 'generic_name = "Tracking2LED"' & key).fetch1('sync_data')

        else:
            sync_data_track = tracking['timestamps']  # [0,1,2,3,4,5]  frames / fs
            sample_rate_frames, img_timestamps = (Tif.SI & (sessions.Session & key)).fetch1(
                'framerate', 'timestamps_sys')

            if np.any(np.diff(img_timestamps) < 0):  # correction for gaps in img-timestamps due to multiple tif files
                sync_data_frames = img_timestamps.copy()
                one_frame = 1 / sample_rate_frames
                gap_ids = np.where(np.diff(sync_data_frames) < 0)[0]
                gap_ids = np.concatenate([gap_ids, [len(sync_data_frames)-1]])
                for seg_start, seg_end in zip(gap_ids[:-1], gap_ids[1:]):
                    sync_data_frames[seg_start+1: seg_end+1] = sync_data_frames[seg_start+1: seg_end+1] + sync_data_frames[seg_start] + one_frame
            else:
                sync_data_frames = img_timestamps

        # Sanity checks
        # 1. Compare length of tracking data and tracking sync data
        # 2. Compare length of spike data and frame sync data 
        # 3. Compare last timestamp frame sync data and tracking sync data

        if len(tracking['x_pos']) != len(sync_data_track):
            raise IndexError('Mismatch between length of sync data and tracking data')
        if len(spikes) != len(sync_data_frames):
            raise IndexError('Mismatch between length of sync data and spiking data')
        if np.abs(sync_data_track[-1] - sync_data_frames[-1]) > np.mean(np.diff(sync_data_frames)):
            raise IndexError('There is more than one frame difference between the end of sync streams')

        # Get seconds to cell to calculate sync sample shift
        seconds_per_plane = 1 / sample_rate_frames / num_planes

        # -> Compensate for the mismatch between timestamp at the beginning of each frame and the time it takes the laser to reach the cell body
        if '2Pmini' in experiment_type:
            seconds_per_line = (Tif.SI & key).fetch1('seconds_per_line')
            seconds_to_cell  = center_plane * seconds_per_plane + center_y * seconds_per_line
            samples_to_cell  = seconds_to_cell * sample_rate_frames # from 'Tif.SI'
            samples_offset   = samples_to_cell + (st_time_offset * sample_rate_frames) # from 'MapParams'
        else:
            raise NotImplementedError('Cell time finding not implemented for experiment type "{}""'.format(experiment_type))


        ######### CREATE SHUFFLING VECTOR ############################################################################

        margin_seconds = shuffle_params['margin_seconds'] 
        break_seconds  = shuffle_params['break_seconds']
        ### TODO: CORRECT THE FOLLOWING LINE ! 
        samples_break = np.ceil(break_seconds / (np.diff(sync_data_frames).mean()/sample_rate_frames)).astype(int)
        sample_from_start = sync_data_frames[0] + (sample_rate_frames * margin_seconds)
        sample_from_end   = sync_data_frames[-1] - (sample_rate_frames * margin_seconds)
        range_start = find_nearest(sync_data_frames, sample_from_start)
        range_end   = find_nearest(sync_data_frames, sample_from_end)

        population_rolling_idxs = range(range_start, range_end, samples_break)
        number_of_shuffles = min(shuffle_params['number_shuffles'], len(population_rolling_idxs))
        shuffling_offsets = random.sample(population_rolling_idxs, number_of_shuffles)


        ######### INSERT INTO MASTER TABLE ##########################################################################    
        key_ = {
            'number_shuffles'  : len(shuffling_offsets), 
            'shuffling_offsets': shuffling_offsets
            }
        self.insert1({**key, **key_})

        # Clean up key
        _ = key.pop('sync_dataset_frames_imaging', None)
        _ = key.pop('sync_name_frames_imaging', None)


        ######### PREPARE SIGNAL AND SHUFFLES ########################################################################

        # Only do this for spikes right now
        signal_data   =  spikes
        signal_idxs   =  np.argwhere(spikes > 0).squeeze()

        # Shuffled array
        shuffled_bvs = []

        ######### SHUFFLE ###########################################################################################

        for shift in tqdm(shuffling_offsets):
            rolled_sync_data_frames = np.roll(sync_data_frames, shift)

            # Look up signal tracking
            signal_dict = calc_signal_dict(rolled_sync_data_frames, sync_data_track,
                                                    tracking, signal_data, signal_idxs, samples_offset, st_params)

            # Get ratemap
            try:
                ratemap_dict = calc_ratemap(occupancy, x_edges, y_edges, signal_dict, ratemap_params)
            except IndexError:  # This happens if only one signal value was found overall
                continue
            
            ratemap_nans =  np.ma.filled(ratemap_dict['ratemap'], fill_value=np.nan).astype(np.float64)
            try:
                if bvs_field_method == 'opexebo':
                    _, fields_map = opexebo.analysis.place_field(
                                        ratemap_nans, min_bins=field_params['min_bins'],
                                        min_mean=ratemap_dict['ratemap'].max()*field_params['fraction_min_mean'],
                                        min_peak=ratemap_dict['ratemap'].max()*field_params['fraction_min_peak'],
                                        init_thresh=field_params['init_thresh'], search_method=field_params['search_method'])

                    fieldmap = fields_map.copy()
                    fieldmap[fieldmap>0] = 1
                elif bvs_field_method == 'bvs':
                    fieldmap, _ = detect_fields(ratemap_nans, minBin=bvs_field_params['min_bin'], \
                                                            std_include=bvs_field_params['std_include'], \
                                                            std_detect=bvs_field_params['std_detect'],\
                                                            show_plots=False, debug=False)
            except: 
                # Whatever, just catch for now ... 
                continue 
                

            if (fieldmap==0).all():
                # Empty fieldmap! 
                continue 

            bvs_params = (BVScoreParams & key).fetch1()
            bvs, _, _, _ = calc_bv_score(fieldmap, r=bvs_params['r_factor'], \
                                        barwidth_max=bvs_params['barwidth_max'], show_plots=False, debug=False)
            shuffled_bvs.append(bvs)


        ######### FILL PART TABLES ################################################################################## 
    
        BVS_Score_dict = {
            'bvs_99'                :  np.nanpercentile(shuffled_bvs, 99),
            'bvs_95'                :  np.nanpercentile(shuffled_bvs, 95),
            'bvs_shuffles'          :  np.array(shuffled_bvs).astype(float)
        }
        self.BVS.insert1({**key, **BVS_Score_dict}, ignore_extra_fields=True)
        