function [fields,allFields] = fieldDetection2D(map,varargin)
% DETECT FIELDS AND CALCULATE BOUNDARY VECTOR SCORE

% basic field detection mechanism written by Oyvind (like he uses in his
% object vector cells paper)
% but added different way of finding the thresholds for field detection
% a(meaning where it will draw the contours) and field inclusion (meaning
% how high the peak firing rate of the field must be for it to be included)


inp = inputParser();
inp.addParameter('median',[])
inp.addParameter('std',[]);
inp.addParameter('minRate',2); % minimum peak rate (Hz) for a field to be included
inp.addParameter('minBin',16); % minimum size (number of pixles) of fields to be included

inp.addParameter('doPlot',true);
inp.addParameter('addMap',true);
inp.addParameter('ax',[]);
inp.addParameter('r',0.5);
inp.parse(varargin{:});
inp.KeepUnmatched = true;
p = inp.Results;

%% Get thresholds for field detection and inclusion 
% thresholds derived from standard deviation and median of the map

if isempty(p.median)
    mapMedian = nanmedian(map(:));
else
    mapMedian = p.median;
end
if isempty(p.std)
    mapStd = nanstd(map(:));
else
    mapStd = p.std;
end

fieldDetectionThresh = mapMedian + mapStd;
fieldInclusionThresh = mapMedian + 2*mapStd;



%% find the fields with countour

% Zeropad to avoid border effects
map(isnan(map)) = 0;
padsize = 3;
padadd = padsize - 1;
mapP = padarray(map,[padsize padsize],0,'both');
mapP(isnan(mapP)) = 0;

% make 1d map
mapP1d = mapP(:);
maxmap = max(mapP1d);

% Get  contourmatrixs
[X,Y] = meshgrid(1:size(mapP,2),1:size(mapP,1));

contourFig = figure('visible', 'on');
[~,h] = contour(X,Y,mapP);

if fieldDetectionThresh < maxmap
    set(h, 'LevelList', linspace(fieldDetectionThresh,maxmap,20));
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
    x = reshape(X,1,numel(X));
    y = reshape(Y,1,numel(Y));
    cm = get(h,'ContourMatrix');
    index_start = 2;
    index_end = cm(2,1)+1;
    IN = inpolygon(x,y,cm(1,index_start:index_end),cm(2,index_start:index_end));
    for i=2:numel(get(h,'LevelList'))
        index_start = index_end + 2;
        index_end = index_start + cm(2,index_start-1) - 1;
        tmp = inpolygon(x,y,cm(1,index_start:index_end),cm(2,index_start:index_end));
        IN = IN | tmp;
    end
    %close(contourFig)
else
    IN = false(1,numel(h.XData));
end

%%
% Define matrix of fields and find connected regions
nap = mapP1d'.*IN;
nap = reshape(nap,size(mapP));
CC = bwconncomp(nap);

% Remove fields that have fewer pixels or lower peak rate than threshold
H = CC.PixelIdxList;
H = H(cellfun(@(x) length(x) >= p.minBin, H));

% wth
H = H(cellfun(@(x) max(nap(x)) >= fieldInclusionThresh, H));

%?
% Remake map with fields all above threshold
fieldMap = zeros(size(mapP1d,1),1);

for i = 1:size(H,2)
    fieldMap(H{1, i}) = mapP1d(H{1, i});
end

fieldMap = reshape(fieldMap,size(mapP));

% Remove zeropad
fieldMap(1:padsize,:) = [];
fieldMap(end-padadd:end,:) = [];
fieldMap(:,1:padsize) = [];
fieldMap(:,end-padadd:end) = [];

% Define connected regions and get centroids
CC = bwconncomp(fieldMap);
S = regionprops(CC,'centroid');


%%  create output

fields = struct('pixelIdxList',{},'incAbThresh',{},'subX',{},'subY',{},'map',{},...
    'centroid',{},'xScore',{},'yScore',{});
for i = 1:CC.NumObjects
    fields(i).pixelIdxList = CC.PixelIdxList{i};
    fields(i).thresh = fieldDetectionThresh;
    
    [sX,sY] = ind2sub(size(map),CC.PixelIdxList{i});
    fields(i).subX = sX;
    fields(i).subY = sY;
    
    tmpFieldMap = zeros(size(fieldMap));
    tmpFieldMap(sX,sY) = fieldMap(sX,sY);
    fields(i).map = tmpFieldMap;
    
    fields(i).incAbThresh = (max(tmpFieldMap)-fieldDetectionThresh)/fieldDetectionThresh;
    fields(i).centroid = S(i).Centroid;
end


if ~p.addMap
    for i = 1:numel(fields)
        fields(i).map = [];
    end
end


allFields = struct('map',[],'numFields',[],...
    'xScore',[],'yScore',[]);
if p.addMap
    allFields.map = fieldMap;
end
allFields.numFields = numel(fields);
[xScore,yScore] = getScores(fieldMap,p.doPlot,p.ax,p.r);
allFields.xScore = xScore;
allFields.yScore = yScore;

allFields.mapMedian = mapMedian;
allFields.mapStd = mapStd;


% Plot map with fields (if doPlot)
if p.doPlot
    if isempty(p.ax)
        ax1 = subplot(2,2,1); 
        ax2 = subplot(2,2,2);
    else
        ax1 = p.ax(1);
        ax2 = p.ax(2);
    end
    imagesc(ax1,map);
    axis(ax1,'square','off','xy');

    imagesc(ax2,fieldMap);
    axis(ax2,'square','off','xy');

    colormap('jet');
end

end


