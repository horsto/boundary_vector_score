### Debora's field detection 
import numpy as np

from skimage import measure, morphology
from scipy import ndimage


# for plotting 
from matplotlib import pyplot as plt
import seaborn as sns


def detect_fields(rmap, std_detect=1, std_include=2, minBin=16, show_plots=False, debug=False):
    ''' 
    Detect fields
    Debora Ledergerber (dlederger@ethz.ch) & Øyvind Arne Hoydal 
    
    Parameters
    ----------
    rmap             :  2-dim np.array
                        (Unprocessed) ratemap. Zero occupancy (nan) are converted to zeros
    std_detect       :  float
                        Number of standard deviations over median for field detection 
    std_include      :  float 
                        Number of standard deviations over median for field inclusion
    minBin           :  integer
                        Minimum no of bins for field inclusion
    show_plots       :  bool 
                        Draw figure? 
    debug            :  bool
                        Show debug messages?
           
    Returns 
    -------
    
    fieldmap         : 2-dim np.array
                       Array in shape of rmap, where detected fields = 1, else 0
    remaining_fields : np.array
                       skimage.measure regionprops output of filtered (remaining) fields
    
    '''
    assert isinstance(rmap,np.ndarray) and len(rmap.shape) == 2, 'Please feed in a 2 dimensional numpy array'
    assert minBin >= 0, 'Parameter "minBin" has to be greater or equal 0'
    assert std_detect >= 1, 'Parameter "std_detect" has to be greater or equal 1'
    assert std_include >= 1, 'Parameter "std_include" has to be greater or equal 1'
    
    mapMedian = np.nanmedian(rmap)
    mapStd    = np.nanstd(rmap)
    # Convert nans to zeros
    ratemap   = np.nan_to_num(rmap)
    # Get thresholds
    fieldDetectionThresh = mapMedian + std_detect * mapStd
    fieldInclusionThresh = mapMedian + std_include * mapStd
    if debug: 
        print('Thresholds\n----------')
        print(f'fieldDetectionThresh: {fieldDetectionThresh}\nfieldInclusionThresh: {fieldInclusionThresh}\nminBin: {minBin}\n')

    ratemap_det = ratemap.copy()
    ratemap_det[ratemap_det<fieldDetectionThresh] = 0
    
    # Convert field map to binary map [0,1]
    thresh_map = ratemap_det.copy()
    thresh_map[thresh_map>0] = 1

    # Find connected components
    labeled_img  = morphology.label(thresh_map, connectivity=2) 

    # Measure detected fields
    regions = measure.regionprops(labeled_img)
    regions = np.array(regions)

    if debug: 
        print('Fields\n------')
        print(f'Detected {len(regions)} fields')

    # Extract maxima and lengths
    maxima_regions = np.array([np.max(ratemap_det[regions[i].coords[:,0], regions[i].coords[:,1]]) for i in range(len(regions))])
    len_regions    = np.array([len(regions[i].coords) for i in range(len(regions))])
    
    if debug: 
        print(f'Maxima of detected regions:   {maxima_regions}')
        print(f'Length of detected regions:   {len_regions}')

    # Apply filtering
    filtered_field_idxs = np.where((len_regions>=minBin) & (maxima_regions>=fieldInclusionThresh))[0]
    remaining_fields    = regions[filtered_field_idxs]

    if debug:
        print(f'{len(remaining_fields)} fields remain\n')
    
    # Create a map from the remaining fields
    fieldmap = np.zeros_like(ratemap_det)
    for field in remaining_fields:
        fieldmap[field.coords[:,0],field.coords[:,1]] = 1
        
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
        ax3.imshow(fieldmap)
        ax3.set_title('Filtered fields')
        sns.despine(left=True,bottom=True)
        
    return fieldmap, remaining_fields


