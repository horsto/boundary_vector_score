function [xScore,yScore] = getScores(map,doPlot,ax,r)
    % Calculate boundary vector score from map
    % r is usually 0.5
    %
    % Returns xScore and yScore (score for each bar
    % orientation)
    
    xDim = size(map,2);
    yDim = size(map,1);
    %r = 0.5;
    
    % map operations
    mapOnes = zeros(size(map));
    mapOnes(map > 1) = 1;
    nMap = numel(map);
    
    % x dim rectangle overlay
    noOvl = nan(5,yDim-5);
    isOvl = nan(5,yDim-5);
    %for i = 1:1 % this would be to vary the size of the mask 
    i = 1; 
    for j = 1:yDim-i
        dummyMap = zeros(size(map));
        dummyMap(j:j+i-1,1:xDim) = 1;
        nDum = numel(find(dummyMap == 1));
        sumMaps = (dummyMap + mapOnes);
        noOvl(i,j) = numel(find(sumMaps == 1))/(nMap-nDum); %
        isOvl(i,j) = numel(find(sumMaps == 2))/nDum; %
    end
    matchScore = isOvl - r*noOvl;
    %end
    [maxJ,idx] = max(matchScore);
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
    noOvl = nan(5,xDim-5);
    isOvl = nan(5,xDim-5);
    for i = 1:1
        for j = 1:xDim-i
            dummyMap = zeros(size(map));
            dummyMap(1:yDim,j:j+i-1) = 1;
            nDum = numel(find(dummyMap == 1));
            sumMaps = (dummyMap + mapOnes);
            noOvl(i,j) = numel(find(sumMaps == 1))/(nMap-nDum); %
            isOvl(i,j) = numel(find(sumMaps == 2))/nDum; %
        end
    end
    matchScore = isOvl - r*noOvl;
    [maxJ,idx] = max(matchScore);
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
        %title('X dimension')
        
        imagesc(ax2,ydMap);
        axis(ax2,'square','off','xy');
        %title('Y dimension')
    end
end

