import datajoint as dj
import numpy as np 


# Load base schema
schema = dj.schema(dj.config['dj_imaging.database'])
schema.spawn_missing_classes()

# Load personal schema 
borderscore_schema = dj.schema('user_horsto_borderscore')
borderscore_schema.spawn_missing_classes()

from bvs.detect_fields import detect_fields
from bvs.bv_score import calc_bv_score


@borderscore_schema 
class BVFieldParams(dj.Lookup):
    definition = """
    # Boundary vector score (BVS) field detection params
    bvfield_params_id        : char(1)        # Parameter set ID, starting with A
    ---       
    std_detect               : float          # Number of standard deviations over median for field detection 
    std_include              : float          # Number of standard deviations over median for field inclusion
    min_bin                  : smallint       # Minimum no of bins for field inclusion
    """  
    contents = [
          ['A', 1., 2., 16]
               ]


@borderscore_schema
class BVField(dj.Computed):
    definition = """
    # Boundary vector score (BVS) fields
    -> Ratemap
    -> BVFieldParams
    ---
    no_fields  = NULL              : smallint        # Total number of detected fields
    fields_map = NULL              : blob@imgstore   # Map of detected fields 
    """
    class Fields(dj.Part):
        definition = """
        # Field detection results
        -> master
        field_no                   : int             # Field number
        ---
        field_coords               : blob@imgstore   # Coordinates of all bins in the firing field
        field_centroid_x           : double          # Field centroid x coordinate
        field_centroid_y           : double          # Field centroid y coordinate
        field_area                 : int             # Area in number of bins
        field_bbox                 : blob@imgstore   # Field bounding box
        """
    @property
    def key_source(self):
        return super().key_source  & 'field_params_id = "A"' 
        # We just want a ratemap + do not care about field detection 
    
    def make(self, key):
        params  =  (BVFieldParams & key).fetch1()
        ratemap_entry = (Ratemap & key).fetch1()
        
        # Process ratemap 
        rm      = np.ma.array(ratemap_entry['ratemap'], mask = ratemap_entry['mask_rm'])
        rm_nans = np.ma.filled(rm, fill_value=np.nan).astype(np.float64)
        
        key['fields_map'], remaining_fields = detect_fields(rm_nans, minBin=16, show_plots=False, debug=False)
        key['no_fields'] = len(remaining_fields)
        
        self.insert1(key)
        
        # Fill "Fields" part table
        for no,field in enumerate(remaining_fields):
            fields_dict = {
                'field_no'         : no,
                'field_coords'     : field.coords,
                'field_centroid_x' : field.centroid[0],
                'field_centroid_y' : field.centroid[1],
                'field_area'       : field.area,
                'field_bbox'       : field.bbox
            }
            self.Fields.insert1({**key,**fields_dict}, ignore_extra_fields=True)



@borderscore_schema 
class BVScoreParams(dj.Lookup):
    definition = """
    # Boundary vector score (BVS) analysis params
    bvscore_params_id        : char(1)        # Parameter set ID, starting with A
    ---       
    r_factor                 : float          # r-factor (weighing contribution of extra fields in diminishing score)
    barwidth_max             : tinyint        # Maximum bar width 
    """  
    contents = [
          ['A', .5, 5]
               ]

@borderscore_schema 
class BVScoreFieldMethod(dj.Lookup):
    definition = """
    # Field detection method for calculating boundary vector score (bvs)
    bv_field_dect_method     : enum('opexebo','bvs')  # Specifies how fields were extracted
    """  
    contents = [
      ['opexebo'],
      ['bvs']
           ]

@borderscore_schema
class BVScore(dj.Computed):
    definition = """
    # Boundary vector score (BVS)
    -> Ratemap
    -> BVField
    -> BVScoreFieldMethod
    ---
    bvs   = NULL               : double         # Boundary vector score (BVS)
    orientation = NULL         : enum('vertical','horizontal')
    """
    class Score_X(dj.Part):
        definition = """
        # Score x (bars spanning X)
        -> master
        ---
        score_x        :  double   # Maximum score for bars spanning x (horizontal bars)
        bar_width      :  tinyint  # Barwidth at maximum 
        ypos           :  smallint # Bar position at maximum (center of bar)
        ypos_rel       :  float    # relative bar position at maximum (center of bar)
        bar_map        :  blob@imgstore   # barMap (streak of ones) at maximum 
        """
    class Score_Y(dj.Part):
        definition = """
        # Score y (bars spanning Y)
        -> master
        ---
        score_y        :  double   # Maximum score for bars spanning y (vertical bars)
        bar_width      :  tinyint  # Barwidth at maximum 
        xpos           :  smallint # Bar position at maximum (center of bar)
        xpos_rel       :  float    # relative bar position at maximum (center of bar)
        bar_map        :  blob@imgstore   # barMap (streak of ones) at maximum 
        """   
    
    @property
    def key_source(self):
        return super().key_source  & 'field_params_id = "A"' 
        # We do not care about other field detection parameters in opexebo right now
        
    def make(self, key):
        # Implement two methods for field detection 
        # 1. opexebo (raw ratemap fieldmaps entry)
        # 2. bvfield 
        
        method = (BVScoreFieldMethod & key).fetch1('bv_field_dect_method')
        if method not in ['opexebo','bvs']:
            raise NotImplementedError(f'Method "{method}" does not have a matching score calculation routine')
        
        try:
            if method == 'opexebo':
                fieldmap = (Ratemap & key).fetch1('fields_map')
                fieldmap = fieldmap.copy()
                fieldmap[fieldmap>0] = 1
            elif method == 'bvs':
                fieldmap = (BVField & key).fetch1('fields_map')
        except: 
            # Whatever, just catch for now ... 
             self.insert1(key)
             return 
             
        if (fieldmap==0).all():
            # Empty fieldmap! 
            self.insert1(key)
            return 

        params = (BVScoreParams & key).fetch1()
        key['bvs'], bvs_x, bvs_y, key['orientation'] = calc_bv_score(fieldmap, r=params['r_factor'], \
                                                barwidth_max=params['barwidth_max'], show_plots=False, debug=False)
        
        self.insert1(key)
        # Part tables
        self.Score_X.insert1({**key,**bvs_x}, ignore_extra_fields=True)
        self.Score_Y.insert1({**key,**bvs_y}, ignore_extra_fields=True)
           