### Debora's borderscore 
import numpy as np

from skimage import measure, morphology
from scipy import ndimage


# for plotting 
from matplotlib import pyplot as plt
import seaborn as sns


def detect_fields(rmap, std_detect=1, std_include=2, minBin=16, show_plots=True):
    ''' 
    Detect fields
    
    Parameters
    ----------
    rmap         :  2-dim np.array
                    (Unprocessed) ratemap. Zero occupancy (nan) are converted to zeros
    std_detect   :  float
                    Number of standard deviations over median for field detection 
    std_include  :  float 
                    Number of standard deviations over median for field inclusion
    minBin       :  integer
                    minimum no of bins for field inclusion
           
    Returns 
    -------
    
    mapOnes      : 2-dim np.array
                   Array in shape of rmap, where detected fields = 1, else 0
    regions      : skimage.measure regionprops output of filtered (remaining) fields
    
    '''
    assert isinstance(rmap,np.ndarray) and len(rmap.shape) == 2, 'Please feed in a 2 dimensional numpy array'
    assert minBin >= 0, 'Parameter "minBin" has to be greater or equal 0'
    assert std_detect >= 1, 'Parameter "std_detect" has to be greater or equal 1'
    assert std_include >= 1, 'Parameter "std_include" has to be greater or equal 1'
    
    mapMedian = np.nanmedian(rmap)
    mapStd    = np.nanstd(rmap)
    # Convert nans to zeros
    ratemap   = np.nan_to_num(rmap)
    
    fieldDetectionThresh = mapMedian + std_detect * mapStd
    fieldInclusionThresh = mapMedian + std_include * mapStd
    
    ratemap_det = ratemap.copy()
    ratemap_det[ratemap_det<fieldDetectionThresh] = 0
    
    thresh_map = ratemap_det.copy()
    thresh_map[thresh_map>0] = 1
    labels = morphology.label(thresh_map, connectivity=2)
    
    # Measure detected fields
    regions = measure.regionprops(labels)
    regions = np.array(regions)

    # Extract maxima and lengths
    maxima_regions = np.array([np.max(ratemap_det[regions[i].coords[:,0], regions[i].coords[:,1]]) for i in range(len(regions))])
    len_regions    = np.array([len(regions[i].coords) for i in range(len(regions))])
    
    # Apply filtering
    filtered_field_idxs = np.where((len_regions>minBin) & (maxima_regions>fieldInclusionThresh))[0]
    remaining_fields    = regions[filtered_field_idxs]
    
    # Create a map from the remaining fields
    mapOnes = np.zeros_like(ratemap_det)
    for no,field in enumerate(remaining_fields):
        mapOnes[field.coords[:,0],field.coords[:,1]] = 1
        
    if show_plots:
        # Generate figure
        figure = plt.figure(figsize=(15,5))
        ax1 = figure.add_subplot(131)
        ax1.imshow(ratemap)
        ax1.set_title('Original ratemap')
        ax2 = figure.add_subplot(132)
        ax2.imshow(ratemap_det)
        ax2.set_title('Detected fields')
        ax3 = figure.add_subplot(133)
        ax3.imshow(mapOnes)
        ax3.set_title('Filtered fields')
        sns.despine(left=True,bottom=True)
        
    return mapOnes, remaining_fields



def calc_borderscore(fieldmap, r=.5):
    ''' 
    Calculate Debora's borderscore
    
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
    borderscore_x = np.nanmax(matchScore)
    
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
    borderscore_y = np.nanmax(matchScore)
    
    return np.max([borderscore_x, borderscore_y]), borderscore_x, borderscore_y
