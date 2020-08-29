function [xScore,yScore] = getScoresPxl(map) 
    % DEPRECATED CODE!!! 
    % USE MAIN FUNCTION: fieldDetection2D (which calls "getScores")
    % TO DETECT FIELDS + CALCULATE BOUNDARY VECTOR SCORE 
        
    % zeropad map
    mapP = padarray(map,[1 1],0,'both');
    yDim = size(map,2);
    xDim = size(mapP,1);  
    
    % neighbors
    xNeighbors = nan(xDim,1);
    yNeighbors = nan(yDim,1);
    
    diagonal1 = nan(xDim,1);
    diagonal2 = nan(xDim,1); 
    
    % find neighbors in each collumn in x dimention
    for k = 1:xDim
        cCol = logical(mapP(k,:)');
        scIdx = find(diff(cCol)); % find state changes
        xNeighbors(k) = sum(scIdx(2:2:end) - scIdx(1:2:end));
    end
    
    % find neighbors in each collumn in y dimention
    for k = 1:yDim
        cCol = logical(mapP(:,k));
        scIdx = find(diff(cCol)); % find state changes
        yNeighbors(k) = sum(scIdx(2:2:end) - scIdx(1:2:end));
    end
    
    % find neighbors in diagonal 1 and 2 
    for k = 1:xDim-1
        % diagonal 1
        cRow = logical([mapP(k,:),0]);
        nRow = logical([0,mapP(k+1,:)]); % next row - shifted to right
        vDiff = logical(abs(diff([cRow;nRow],1))); % find the difference between the rows
        diagonal1(k) = numel(find(~vDiff & cRow)); % find where there is no difference but there are values in current row
        
        % dagonal 2
        cRow = logical([0,mapP(k,:)]);
        nRow = logical([mapP(k+1,:),0]); % next row - shifted to right
        vDiff = logical(abs(diff([cRow;nRow],1))); % find the difference between the rows
        diagonal2(k) = numel(find(~vDiff & cRow)); % find where there is no difference but there are values in current row
    end
    
    
    xScore = abs(nansum(xNeighbors) - max(nansum(diagonal1),nansum(diagonal2)))/nansum(xNeighbors); 
    yScore = abs(nansum(yNeighbors) - max(nansum(diagonal1),nansum(diagonal2)))/nansum(yNeighbors);

end
