### Debora's boundary vector score 
import numpy as np


def calc_bv_score(fieldmap, r=.5):
    ''' 
    Calculate boundary vector (bv) score
    Debora Ledergerber, dlederger@ethz.ch
    
    Parameters
    ----------
    fieldmap : 2-dim np.array
               (Filtered) field map where in-field == 1, else 0
               
    r        : float (default = 0.5)
               r-factor (weighing contribution of extra fields in diminishing score)
               If 0, other fields do not contribute to score 
               
    TODO: Experiment with i (the width over the overlay bar)
    TODO: Extract position of bar maximum / distance from wall 
    '''
    assert isinstance(fieldmap,np.ndarray) and len(fieldmap.shape) == 2, 'Please feed in a 2 dimensional numpy array'
    assert 0 <= r <= 1, 'Parameter "r" has to be a float in between 0 and 1'
    
    xDim = fieldmap.shape[1]
    yDim = fieldmap.shape[0]
    nMap = np.size(fieldmap)
    
    nan_map = np.zeros_like(fieldmap)
    nan_map[nan_map==0] = np.nan

    # ... in x
    # Create "result bank"
    noOvl = nan_map.copy()
    isOvl = nan_map.copy()

    i = 1
    for j in range(yDim-i):
        # Create a horizontal bar
        dummyMap = np.zeros_like(fieldmap)
        dummyMap[j:j+i,:xDim] = 1
        # ... and count how many ones you created
        nDum    = len(dummyMap[dummyMap==1])

        # Add bar to original map
        sumMaps = (dummyMap + fieldmap)
        noOvl[i,j] = len(sumMaps[sumMaps==1])/(nMap-nDum)
        isOvl[i,j] = len(sumMaps[sumMaps==2])/nDum
    # Score in x
    matchScore = isOvl - r*noOvl
    bvs_x = np.nanmax(matchScore)
    
    # ... in y
    # Create "result bank"
    noOvl = nan_map.copy()
    isOvl = nan_map.copy()
    
    i = 1
    for j in range(xDim-i):
        # Create a vertical bar
        dummyMap = np.zeros_like(fieldmap)
        dummyMap[:yDim, j:j+i] = 1
        # ... and count how many ones you created
        nDum    = len(dummyMap[dummyMap==1])

        # Add bar to original map
        sumMaps = (dummyMap + fieldmap)
        noOvl[i,j] = len(sumMaps[sumMaps==1])/(nMap-nDum)
        isOvl[i,j] = len(sumMaps[sumMaps==2])/nDum
        
    # Score in y
    matchScore = isOvl - r*noOvl
    bvs_y = np.nanmax(matchScore)
    
    return np.max([bvs_x, bvs_y]), bvs_x, bvs_y
