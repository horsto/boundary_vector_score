## HELPERS 
# copied from https://github.com/kavli-ntnu/dj-moser-imaging

import numpy as np 
import math
from astropy.convolution import convolve, Gaussian2DKernel


def find_nearest(array, value):
    ''' Find nearest element to "value" in "array" and return index of that element '''
    idx = np.searchsorted(array, value, side='left')
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx
    
    
def calc_signal_dict(signal_sync, tracking_sync, tracking_data, signal_data, signal_idxs, samples_offset, params):
    '''
    Look up tracking signal for every valid calcium event/spike based on sync data.
    
    Parameters
    ----------
    signal_sync : np.array
        Sync pulse number for every event in 'signal_data'
    tracking_sync : np.array
        Sync pulse number for every event in 'tracking_data'        
    tracking_data : dict
        Fetched Tracking.OpenField() entry with:
            - x_pos
            - y_pos
            - head_angle
            - speed
    signal_data : np.array
        Calcium or spikes signal of length 'signal_sync'
    signal_idxs : np.array 
        Indices of signal to keep (filter)
    samples_offset : int
        Number of samples to add to 'signal_sync' samples before 
        lookup of tracking data. This shifts signal in time compared 
        to tracking data and is used (1) to compensate for laser fly time 
        in FOV and (2) for shuffling purposes 
    params: dict
        - speed_cutoff_low: lower speed cutoff in unit 'speed' has in 'tracking_data'
        - speed_cutoff_low: upper speed cutoff in unit 'speed' has in 'tracking_data'
       
    Returns
    -------
    signal_dict : dict
        x_pos_signal
        y_pos_signal
        head_angle_signal
        speed_signal
        signal (amplitudes)
    '''
    
    tracking_keys = {'x_pos', 'y_pos', 'head_angle', 'speed'}
    signal_dict = {}
    
    tracking_indices_filtered = []  # Speed filtered
    signal_indices_filtered   = []  # Spikes / Fluorescence

    signal_sync          = signal_sync.astype(float)
    signal_sync          += samples_offset  # Shift in time
    signal_sync          = signal_sync[signal_idxs]
    
    for idx, sync_pulse in zip(signal_idxs, signal_sync): # signal_idxs either filtered or whole series (spikes vs. delta f/f)
        tracking_idx = find_nearest(tracking_sync, sync_pulse)
        if (tracking_data['speed'][tracking_idx] > params['speed_cutoff_low']) and (tracking_data['speed'][tracking_idx] < params['speed_cutoff_high']):
            tracking_indices_filtered.append(tracking_idx)
            signal_indices_filtered.append(idx)
            
    # Build signal tracking dictionary
    for pos_key in tracking_keys:
        signal_dict[pos_key + '_signal'] = tracking_data[pos_key][tracking_indices_filtered]

    signal_dict['signal'] = signal_data[signal_indices_filtered].squeeze()
    if not signal_dict["signal"].ndim == 1:
        # Somewhat convoluted logic: I'm not sure what circumstances call for the `.squeeze()`,
        # but where there is only a single value, that compresses it to a 0d array
        # Datajoint later converts that to a non-array floating point number during the 
        # conversion to/from a blob. Therefore, enforce that the array still has at least 1 dimension
        signal_dict["signal"] = np.expand_dims(signal_dict["signal"], 0)
    return signal_dict


def calc_ratemap(occupancy, x_edges, y_edges, signaltracking, params):
    '''
    Calculate ratemap
    Parameters
    ----------
    occupancy : masked np.array
        Smoothed occupancy. Masked where occupancy low
    x_edges : np.array
        Bin edges in x 
    y_edges : np.array
        Bin edges in y
    signaltracking : dict
        SignalTracking table entry
    params : dict
        MapParams table entry
    
    Returns
    -------
    ratemap_dict : dict
        - binned_raw : np.array: Binned raw (unsmoothed) signal
        - ratemap_raw: np masked array: Unsmoothed ratemap (mask where occupancy low)
        - ratemap    : np masked array: Smoothed ratemap (mask where occupancy low)
        - bin_max    : tuple   : (x,y) coordinate of bin with maximum signal
        - max        : float : Max of signal 
        
    '''
    ratemap_dict = {}
    
    binned_signal = np.zeros_like(occupancy.data)
    # Add one at end to not miss signal at borders
    x_edges[-1] += 1
    y_edges[-1] += 1

    # Look up signal per bin
    for no_x in range(len(x_edges)-1):
        for no_y in range(len(y_edges)-1):
            boolean_x = (signaltracking['x_pos_signal'] >= x_edges[no_x]) & (signaltracking['x_pos_signal'] < x_edges[no_x+1])
            boolean_y = (signaltracking['y_pos_signal'] >= y_edges[no_y]) & (signaltracking['y_pos_signal'] < y_edges[no_y+1])
            extracted_signal = signaltracking['signal'][boolean_x & boolean_y]
            binned_signal[no_y, no_x] = np.nansum(extracted_signal)

    ratemap_dict['binned_raw'] = binned_signal
    binned_signal = np.ma.masked_where(occupancy.mask, binned_signal)  # Masking. This step is probably unnecessary
    ratemap_dict['ratemap_raw'] = binned_signal / occupancy
    
    # Instead of smoothing the raw binned spikes, substitute those values that are masked in
    # occupancy map with nans.
    # Then use astropy.convolve to smooth padded version of the spikemap 
        
    binned_signal[occupancy.mask] = np.nan
    kernel = Gaussian2DKernel(x_stddev=params['sigma_signal'])

    pad_width = int(5*params['sigma_signal'])
    binned_signal_padded = np.pad(binned_signal, pad_width=pad_width, mode='symmetric')  # as in BNT
    binned_signal_smoothed = convolve(binned_signal_padded, kernel, boundary='extend')[pad_width:-pad_width, pad_width:-pad_width]
    binned_signal_smoothed = np.ma.masked_where(occupancy.mask, binned_signal_smoothed)  # Masking. This step is probably unnecessary
    masked_ratemap = binned_signal_smoothed / occupancy

    ratemap_dict['ratemap']       = masked_ratemap
    ratemap_dict['bin_max']       = np.unravel_index(masked_ratemap.argmax(), masked_ratemap.shape)
    ratemap_dict['max']           = np.max(masked_ratemap)
    
    return ratemap_dict
