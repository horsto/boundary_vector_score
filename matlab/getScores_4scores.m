function [xScore,yScore] = getScores_4scores(map)
    % DEPRECATED CODE!!! 
    % USE MAIN FUNCTION: fieldDetection2D (which calls "getScores")
    % TO DETECT FIELDS + CALCULATE BOUNDARY VECTOR SCORE 

    % thresholds
    thrP = 0.8;
    %thrR = 0.1; % divided by extent
    thrE = 0.5;
    thrN = 0.7;
    thrS = 0.7;

    % xDim
    xDim = size(map,2);
    xParallel = nan(xDim,1);
    xExtent = nan(xDim,1);
    % yDim
    yDim = size(map,1);
    yParallel = nan(yDim,1);
    yExtent = nan(yDim,1);
    
    
    % find xParallel (how parallel the field is to the x wall)
    for k = 2:xDim-1
        cCol = map(:,k);
        props = regionprops(true(size(cCol)),cCol,'WeightedCentroid');
        xParallel(k) = props.WeightedCentroid(2);
        tmpExtent = find(logical([0;cCol;0]),1):find(logical([0;cCol;0]),1,'last');
%       tmpExtent = diff(find(diff(logical([0;cCol;0]))~=0));
        if ~isempty(tmpExtent)
            xExtent(k) = numel(tmpExtent)/yDim;
%           xExtent(k) = max(tmpExtent(1:2:end))/yDim;
        end
    end
    
    % find b (y achsenabschnitt und steigung) by fitting a lin regression
    fitX = [ones(1,xDim);1:xDim]';
    vNan = isnan(xParallel);
    fitX(vNan,:) = [];
    xParallel(vNan,:) = [];
    b = fitX\xParallel;
    
    % find R2 (goodness of fit)
    xCalc = fitX*b;
    v1 = [fitX(1,2),b(1)+fitX(1,2)*b(2),0];
    v2 = [fitX(end,2),b(1)+fitX(end,2)*b(2),0];
    for i = 1:numel(xParallel)
        pt = [fitX(i,2),xParallel(i,:),0];
        dm(i) = point_to_line_matlab(pt, v1, v2);
        d(i) = point_to_line(pt, b(1), b(2));
    end
%     Rsq2 = 1 - sum((xParallel - xCalc).^2)/sum((xParallel - mean(xParallel)).^2);
    sprd = 1-median(d)/3.2;
    sprd(sprd<0) = 0;
    % find narrowness of field
    xExtent(~(xExtent>0)) = [];
    
    xScore.closeToWall = abs((b(1)-(yDim/2))/(yDim/2)); % closeness to wall score: from 0-1
    xScore.ptw = abs(1-abs(b(2))); % parralell to wall score: from 0-1
%   xScore.rsq = Rsq2;
    xScore.ext = numel(xParallel)/xDim; % extent of coverage: from 0-1
    xScore.nrw = 1-median(xExtent);
    xScore.sprd = sprd;
    
    if xScore.ptw>thrP && xScore.sprd>thrS && xScore.ext>thrE && xScore.nrw>thrN
        xScore.bc = true;
    else
        xScore.bc = false;
    end
    gcf;
    hold on
    subplot(2,2,3);
    plot(fitX(:,2),xParallel,'k*');
    hold on
    plot(fitX(:,2),xCalc)
    ax = gca;
    ax.XLim = [0,32];
    ax.YLim = [0,32];
    title('X dimension')
  
    

    for k = 2:yDim-1
        cCol = map(k,:)';
        props = regionprops(true(size(cCol)),cCol,'WeightedCentroid');
        yParallel(k) = props.WeightedCentroid(2);
        tmpExtent = find(logical([0;cCol;0]),1):find(logical([0;cCol;0]),1,'last');
%         tmpExtent = diff(find(diff(logical([0;cCol;0]))~=0));
        if ~isempty(tmpExtent)
            yExtent(k) = numel(tmpExtent)/xDim;
%             yExtent(k) = max(tmpExtent(1:2:end))/xDim;
        end
        
    end
    
    % find b (y achsenabschnitt und steigung) by fitting a lin regression
    fitX = [ones(1,yDim);1:yDim]';
    vNan = isnan(yParallel);
    fitX(vNan,:) = [];
    yParallel(vNan,:) = [];
    b = fitX\yParallel;

    % find R2 (goodness of fit)
    yCalc = fitX*b;
    v1 = [fitX(1,2),b(1)+fitX(1,2)*b(2),0];
    v2 = [fitX(end,2),b(1)+fitX(end,2)*b(2),0];
    d = [];
    for i = 1:numel(yParallel)
        pt = [fitX(i,2),yParallel(i,:),0];
        dm(i) = point_to_line_matlab(pt, v1, v2);
        d(i) = point_to_line(pt, b(1), b(2));
    end
%     Rsq2 = 1 - sum((yParallel - yCalc).^2)/sum((yParallel - mean(yParallel)).^2);
    sprd = 1-mean(d)/3.2;
    sprd(sprd<0) = 0;
    ySprd = 1-(mean(abs(yParallel-yCalc))/5);
    
     % find narrowness of field
    yExtent(~(yExtent>0)) = [];

    yScore.closeToWall = abs((b(1)-(xDim/2))/(xDim/2)); % closeness to wall score: from 0-1
    yScore.ptw = abs(1-abs(b(2))); % parralell to wall score: from 0-1
%     yScore.rsq = Rsq2;
    yScore.ext = numel(yParallel)/yDim; % extent of coverage: from 0-1
    yScore.nrw = 1-mean(yExtent);
    yScore.sprd = sprd;
    if yScore.ptw>thrP && yScore.sprd>thrS && yScore.ext>thrE && yScore.nrw>thrN
        yScore.bc = true;
    else
        yScore.bc = false;
    end
    
    
    subplot(2,2,4);
    hold on
    plot(fitX(:,2),yParallel,'k*');
    hold on
    plot(fitX(:,2),yCalc)
    ax = gca;
    ax.XLim = [0,32];
    ax.YLim = [0,32];
    title('Y dimension')
end



function d = point_to_line_matlab(pt, v1, v2)
  a = v1 - v2;
  b = pt - v2;
  d = norm(cross(a,b)) / norm(a);
end

function d = point_to_line(pt, a, b)
  d = abs(pt(2)-(a+b*pt(1)))*cos(atan(b));
end



