function [xScore,yScore] = getScoresNT(map,doPlot,ax)
    % DEPRECATED CODE!!! 
    % USE MAIN FUNCTION: fieldDetection2D (which calls "getScores")
    % TO DETECT FIELDS + CALCULATE BOUNDARY VECTOR SCORE 

    map = map/max(max(map));
    
    xDim = size(map,2);
    yDim = size(map,1);
    r = 1;
    
    % map operations
    mapOnes = zeros(size(map));
    mapOnes(map > 1) = 1;
    nMap = numel(map);
    
    % x dim rectangle overlay
    scores = nan(5,yDim-5);
    for i = 1:5
        for j = 1:yDim-i
            dummyMap = zeros(size(map));
            dummyMap(j:j+i-1,1:xDim) = true;
            
            inBarIdx = find(dummyMap);
            inBarMean = mean(map(inBarIdx));
            
            outBarIdx = find(~dummyMap);
            outBarMean = mean(map(outBarIdx));

            scores(i,j) = inBarMean - r*outBarMean;
        end
    end
    [maxJ,idx] = max(scores);
    [maxMS,imJ] = max(maxJ);
    imI = idx(imJ);
    
    xScore.sc = maxMS;
    xScore.recWidth = imI;
    xScore.pxlIdx = imJ;
    xdMap = zeros(size(map));
    xdMap(imJ:imJ+imI-1,:) = 1;
    xScore.dMap = xdMap;
    
    [row,~] = find(xdMap,1);
    xScore.dist = row/yDim;
    
    % y dim rectangle overlay
    scores = nan(5,yDim-5);
    for i = 1:5
        for j = 1:xDim-i
            dummyMap = zeros(size(map));
            dummyMap(1:yDim,j:j+i-1) = 1;
            
            inBarIdx = find(dummyMap);
            inBarMean = mean(map(inBarIdx));
            
            outBarIdx = find(~dummyMap);
            outBarMean = mean(map(outBarIdx));

            scores(i,j) = inBarMean - r*outBarMean;
        end
    end
    [maxJ,idx] = max(scores);
    [maxMS,imJ] = max(maxJ);
    imI = idx(imJ);
    
    yScore.sc = maxMS;
    yScore.recWidth = imI;
    yScore.pxlIdx = imJ;
    ydMap = zeros(size(map));
    ydMap(:,imJ:imJ+imI-1) = 1;
    yScore.dMap = ydMap;
    
    [~,col] = find(ydMap,1);
    yScore.dist = col/xDim;
    
            
    if doPlot
        if isempty(ax)
            ax1 = subplot(2,2,3);
            ax2 = subplot(2,2,4);
        else
            ax1 = ax(3);
            ax2 = ax(4);
        end
        imagesc(ax1,xdMap);
        axis(ax1,'square','off','xy');
%         title('X dimension')
        
        imagesc(ax2,ydMap);
        axis(ax2,'square','off','xy');
%         title('Y dimension')
    end
end
