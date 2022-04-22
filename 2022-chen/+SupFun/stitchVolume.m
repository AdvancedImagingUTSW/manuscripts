function [stitchedVolume,volCell]=stitchVolume(volCell,tformCell)

sizeCell=cellfun(@(v) size(v),volCell,'unif',0);
assert(all([(cellfun(@(s) all(sizeCell{1}==s),sizeCell))]));

% Compute the new output limits
for vIdx=1:length(volCell) 
 [xLimitsOut(vIdx,:),yLimitsOut(vIdx,:),zLimitsOut(vIdx,:)] = outputLimits(tformCell{vIdx},[1,size(volCell{vIdx},2)],...
     [1 , size(volCell{vIdx},1)],[1 size(volCell{vIdx},3)]);
end

volSize=sizeCell{1};
% Estimate The new size of the image
xMin= min([1; xLimitsOut(:)]);
yMin= min([1; yLimitsOut(:)]);
zMin= min([1; zLimitsOut(:)]);
xMax= max([volSize(2); xLimitsOut(:)]);
yMax= max([volSize(1); yLimitsOut(:)]);
zMax= max([volSize(3); zLimitsOut(:)]);
width  = round(xMax-xMin);
height = round(yMax-yMin);
depth  = round(zMax-zMin);
xLimits = [xMin xMax];
yLimits = [yMin yMax];
zLimits = [zMin zMax];

% Define a outputRef and apply transforms
outputRef=imref3d([height width depth], xLimits,yLimits,zLimits);
%vol2=imwarp(vol2,invert(registrationInfo.tform),'OutputView',outputRef);
for vIdx=1:length(volCell)
    volCell{vIdx}=imwarp(volCell{vIdx},tformCell{vIdx},'OutputView',outputRef);
end 

%% Perform stitching
% Collect median intensities for each plane
% use erosion to mask the plane border.
% collect the number of plane imaged on each pixel in <volsMask>
volsMask=zeros(size(volCell{1}));
% Structure element for erosion
[x,y,z] = ndgrid(-3:3);
se = strel(sqrt(x.^2 + y.^2 + z.^2) <=3);

medianIntensities=zeros(1,length(volCell));
warpedMask=cell(1,length(volCell));
for vIdx=1:length(volCell)
    warpedMask{vIdx}=volCell{vIdx}>0;
    warpedMask{vIdx}=imerode(warpedMask{vIdx},se);
    volsMask(warpedMask{vIdx})=volsMask(warpedMask{vIdx})+1; 
    medianIntensities(vIdx)=median(volCell{vIdx}(warpedMask{vIdx}));
end

% scale the pixel to the middle plane intensity
% average the planes that are present on the same pixels
stitchedVolume=zeros(size(volCell{1}));
for vIdx=1:length(volCell)
    weights=medianIntensities(2)./(medianIntensities(vIdx)*max(1,volsMask));
    weights(warpedMask{vIdx}==0)=0; 
    stitchedVolume=stitchedVolume+weights.*double(volCell{vIdx});
end
