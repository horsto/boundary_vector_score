### Debora's boundary vector score 
import numpy as np

# for plotting 
from matplotlib import pyplot as plt
import seaborn as sns

def calc_bv_score(fieldmap, r=.5, barwidth_max=1, show_plots=False, debug=False):
    ''' 
    Calculate boundary vector (bv) score
    Debora Ledergerber, dlederger@ethz.ch
    Python integration by Horst Obenhaus (horst.obenhaus@ntnu.no)
    
    Parameters
    ----------
    fieldmap     : 2-dim np.array
                   (Filtered) field map where in-field == 1, else 0
               
    r            : float (default = 0.5)
                   r-factor (weighing contribution of extra fields in diminishing score)
                   If 0, other fields do not contribute to score 
    barwidth_max : int 
                   Maximum bar width. Defaults to 1 (Only the most narrow fields get high scores) 
               
    Returns 
    -------
    max_score       : float
                      Maximum of x (horizontal bars) and y (vertical bars) score 
    results_x       : dictionary
                      'score_x'   : Maximum score for bars spanning x (horizontal bars)
                      'bar_width' : barwidth at maximum 
                      'ypos'      : bar position at maximum (center of bar)
                      'ypos_rel'  : relative bar position at maximum (center of bar)
                      'bar_map'   : barMap (streak of ones) at maximum 
    results_y       : dictionary
                      'score_y'   : Maximum score for bars spanning y (vertical bars)
                      'bar_width' : barwidth at maximum 
                      'xpos'      : bar position at maximum (center of bar)
                      'xpos_rel'  : relative bar position at maximum (center of bar)
                      'bar_map'   : barMap (streak of ones) at maximum 
    max_orientation : string or None
                      Maximum orientation, either 'horizontal' or 'vertical' 
                      Returns None if exactly the same in both directions
    show_plots      : bool 
                      Draw figure? 
    debug           : bool
                      Show debug messages?
   
    '''   
    assert isinstance(fieldmap,np.ndarray) and len(fieldmap.shape) == 2, 'Please feed in a 2 dimensional numpy array for "fieldmap"'
    assert (np.max(fieldmap) == 1) and (np.min(fieldmap) == 0), 'Please make sure "fieldmap" is a 2D array where field coordinates = 1 and rest = 0' 
    assert 0 <= r <= 1, 'Parameter "r" has to be a float in between 0 and 1'
    assert isinstance(barwidth_max, int) and (barwidth_max > 0), 'Parameter "barwidth_max" has to be a positive integer value'
    
    xDim = fieldmap.shape[1]
    yDim = fieldmap.shape[0]
    nMap = np.size(fieldmap)
    barwidth_max += 1

    if debug: 
        print(f'Barwidth ranges:  {np.arange(1, barwidth_max)}')
        print(f'r value: {r}')
    # Logic of score calculation:

    # Create maps that contains streaks of ones in x and y 
    # Loop those separately over the field map and add them to the map 
    # Calculate the percentage of overlap of any field with the bar (those have 
    # the value two) and compare with the rest (those have value one)
    # The value "r" is a factor that weighs the contribution of fields not 
    # overlapping with the bar in diminishing the score. Its default value is 0.5.

    # Result that maximise this calculation are saved in 
    # "results_x" and "results_y"

    #### Bars spanning X ############################################################################
    if debug:
        print('Looping over horizontal bars ... ')

    score_max = -1
    results_x = {}
    for barwidth in np.arange(1, barwidth_max):
        for yPos in range(yDim - barwidth + 1):
            
            barMap = np.zeros_like(fieldmap)
            barMap[yPos:yPos+barwidth,:] = 1
            nBar = barwidth * xDim # Number of ones in bar
            
            # Add bar to original map
            sumMaps = (fieldmap + barMap)
            perc_ones = len(sumMaps[sumMaps==1])/(nMap-nBar)
            perc_twos = len(sumMaps[sumMaps==2])/nBar
            score = perc_twos - r*perc_ones
            
            if score > score_max:
                if debug: 
                    print(f'Score increased to {score: .3f} | yPos: {yPos} (barwidth: {barwidth})')

                score_max = score
                results_x['score_x']    = score
                results_x['bar_width']  = barwidth
                results_x['ypos']       = yPos + (barwidth/2)
                results_x['ypos_rel']   = results_x['ypos']/yDim
                results_x['bar_map']    = barMap   

    #### Bars spanning Y ############################################################################
    if debug:
        print('\nLooping over vertical bars ... ')

    score_max = -1
    results_y = {}
    for barwidth in np.arange(1, barwidth_max):
        for xPos in range(xDim - barwidth + 1):
            barMap = np.zeros_like(fieldmap)
            barMap[:, xPos:xPos+barwidth] = 1
            nBar = barwidth * yDim # Number of ones in bar
            
            # Add bar to original map
            sumMaps = (fieldmap + barMap)
            perc_ones = len(sumMaps[sumMaps==1])/(nMap-nBar)
            perc_twos = len(sumMaps[sumMaps==2])/nBar
            score = perc_twos - r*perc_ones
            
            if score > score_max:
                if debug: 
                    print(f'Score increased to {score: .3f} | xPos: {xPos} (barwidth: {barwidth})')

                score_max = score
                results_y['score_y']   = score
                results_y['bar_width'] = barwidth
                results_y['xpos']      = xPos + (barwidth/2)
                results_y['xpos_rel']  = results_y['xpos']/xDim
                results_y['bar_map']   = barMap   


    if show_plots:
        # Draw figure of barMaps at maximum score for horizontal and vertical 
        
        sns.set(style='white', font_scale=1.2)
        plt.rcParams['xtick.major.size'] = 5
        plt.rcParams['xtick.major.width'] = 1
        plt.rcParams['ytick.major.size'] = 5
        plt.rcParams['ytick.major.width'] = 1
        plt.rcParams['xtick.bottom'] = True
        plt.rcParams['ytick.left'] = True

        figure = plt.figure(figsize=(15,5))
        ax = figure.add_subplot(131)
        ax.imshow(fieldmap)
        ax.set_title('Field map')
        ax = figure.add_subplot(132)
        ax.imshow(results_x['bar_map'])
        ax.set_title('Horizontal bar max: {:.2f}'.format(results_x['score_x']))
        ax = figure.add_subplot(133)
        ax.imshow(results_y['bar_map'])
        ax.set_title('Vertical bar max: {:.2f}'.format(results_y['score_y']))
        sns.despine(left=True, bottom=True)

    max_score = np.max([results_x['score_x'], results_y['score_y']])
    
    if results_x['score_x'] > results_y['score_y']:
        max_orientation = 'horizontal'
    elif results_x['score_x'] < results_y['score_y']:
        max_orientation = 'vertical'
    else: 
        max_orientation = None

    if debug:
        print(f'\nFinal boundary vector score: {max_score} ({max_orientation})')
    

    return max_score, results_x, results_y, max_orientation
