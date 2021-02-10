%% Init 
if(isempty(which('MovieData')))
    addpath(genpath('~/repo/utsw-ssh/'));
else
    disp('Code is already loaded');
end

clc;
clear; close all


%% Detect objectes on two projections angle 
data_source = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPM/projection/MV3fixed/GEMS/201002/Cell10/averagedframes/';
proj=double(imread([data_source filesep '/ch1/0000.tif']));
cellProj{1}=double(imread([data_source filesep '/ch1/0000.tif']));
cellProj{2}=double(imread([data_source filesep '/ch1/0019.tif']));
objectScale=1.8;
detectSensitivity=0.01;

cellDetect=cell(size(cellProj));
for pIdx=1:numel(cellProj)
    pts=pointSourceDetection(cellProj{pIdx},objectScale,'Alpha',detectSensitivity);
    
    figure();
    imshow(mat2gray(cellProj{pIdx}));
    hold on;
    scatter(pts.x,pts.y);
    drawnow; 
    hold off
    cellDetect{pIdx}=[pts.x' pts.y'];
end

%% Match objects to the first frame 
close all
matchingSearchRadius=100;
xGating=3; % horizontal gating used for matching 
cellMatchIdx=cell(1,numel(cellDetect)-1);
for pIdx=2:numel(cellProj)
    X1=cellDetect{1};
    X2=cellDetect{pIdx};
    D = createSparseDistanceMatrix(X1,X2, matchingSearchRadius);
    % only keep the vertical distance 
    [i,j,v] = find(D);
    gatedDistanceIdx=(abs(X1(i,1)-X2(j,1))>xGating);
    gatedDistance=sparse(i(~gatedDistanceIdx),j(~gatedDistanceIdx),...
        v(~gatedDistanceIdx),size(D,1),size(D,2));    
    
    [link12, ~] = lap(gatedDistance, [], [], 1);
    
    n1 = size(X1,1);
    n2 = size(X2,1);
    link12 = link12(1:n1);
    matchIdx = link12<=n2;
    idx1 = find(matchIdx);
    idx2 = double(link12(matchIdx));

    cmap=prism(numel(idx1));
       
    figure();
    img=sc(cellProj{1},'summer')+sc(cellProj{pIdx},'autumn');
    sc(img);
    hold on;
    scatter(cellDetect{1}(:,1),cellDetect{1}(:,2),50,[0 0 0]);
    scatter(cellDetect{1}(idx1,1),cellDetect{1}(idx1,2),50,cmap);
    scatter(cellDetect{pIdx}(:,1),cellDetect{pIdx}(:,2),50,[0.2 0.2 0.2]);
    scatter(cellDetect{pIdx}(idx2,1),cellDetect{pIdx}(idx2,2),50,cmap);
    for l=1:numel(idx1)
    plot([cellDetect{1}(idx1(l),1) cellDetect{pIdx}(idx2(l),1)], ...
         [cellDetect{1}(idx1(l),2) cellDetect{pIdx}(idx2(l),2)]);
    end
    hold off;
    title('Naive LAP pairing (non-directional)')

    cellMatchIdx=[idx1 idx2];
end

%% Z estimation
cell3DPos=cell(1,numel(cellMatchIdx));
gain=1:0.01:1.8;  %normalized shear factor (normalized such that 1 equals properly deskewed view)
slice=20; % frame you want to use together with first frame
alpha_2=atand(gain(slice)*2*sind(60)-2*sind(60)); %shear angle
clear h

for pIdx=2:numel(cellProj)
    posRef=cellDetect{1}(idx1,:);
    posMatched=cellDetect{pIdx}(idx2,:);   
    delta_y=posMatched(:,2) - posRef(:,2); 
    Z=delta_y/(tand(alpha_2))+200; % Still do not understand where the 200 comes from
    cell3DPos{pIdx-1}=[posRef Z];
end

XYZ=cell3DPos{1};

%% Crop the huge ground truth
cropFilePath='/tmp/liveProj_crop.tif';
cropOverwrite=false;

if(~exist(cropFilePath,'file')||cropOverwrite)
    % fileGT='/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPM/projection/MV3fixed/GEMS/201002/Cell7/1_CH00_000000_deskewed.tif';
    % stackGTOrig=stackRead(fileGT);
    stackGT=stackRead('/tmp/test.tif');  %% There is bug with the original tiff format, this is resaved with gmic
    % reproducing the mirror and rotationg that was in the original microscope code
    stackGT=flip(stackGT,2);
    stackGT=imrotate3(stackGT,90,[0 0 1]);
    figure();
    imshow(squeeze(max(stackGT,[],3)),[],'Border','tight');
    drawnow;
    % cropTruth reusing Kevin nice trick
    close all
    vol=stackGT;
    maxXY = squeeze(max(vol,[],3));
    figure
    imshow(maxXY,[],'Border','tight');
    xyIndices = ceil(getrect);
    
    vol2 = vol(xyIndices(2):xyIndices(2)+xyIndices(4),xyIndices(1):xyIndices(1)+xyIndices(3),:);
    maxYZ = squeeze(max(vol2,[],1)); imshow(maxYZ,[],'Border','tight');
    zIndices = ceil(getrect);
    
    roiIdx = nan(6,1);
    roiIdx(1) = max(1,xyIndices(1));
    roiIdx(2) = min(size(vol,2),xyIndices(1)+xyIndices(3)-1);
    roiIdx(3) = max(1,xyIndices(2));
    roiIdx(4) = min(size(vol,1),xyIndices(2)+xyIndices(4)-1);
    roiIdx(5) = max(1,zIndices(1));
    roiIdx(6) = min(size(vol,3),zIndices(1)+zIndices(3)-1);
    cropTruth=vol(roiIdx(3):roiIdx(4),roiIdx(1):roiIdx(2),roiIdx(5):roiIdx(6));
    figure
    imshow(maxXY,[],'Border','tight');
    stackWrite(cropTruth,cropFilePath);
else
    cropTruth=stackRead(cropFilePath);
end

figure
imshow(squeeze(max(cropTruth,[],3)),[],'Border','tight'); 

%% Detect object on truth 
pstruct=pointSourceDetection3D(cropTruth,1.8,'Alpha',0.0000001);
close all;
figure();
imshow(squeeze(max(cropTruth,[],3)),[],'Border','tight'); 
hold on;
scatter(pstruct.x,pstruct.y);
hold off; 
drawnow;

figure();
imshow(squeeze(max(cropTruth,[],2)),[],'Border','tight'); 
hold on;
scatter(pstruct.z,pstruct.y);
hold off; 
drawnow;

%% Ground truth and measured are shifted (this is the case even without crop)
XYZMeasured=[pstruct.x',pstruct.y',pstruct.z'];

close all
figure();
subplot(2,2,1);
hold on;
scatter(XYZMeasured(:,1),XYZMeasured(:,2));
scatter(XYZ(:,1),XYZ(:,2));
drawnow;
hold off; 
subplot(2,2,2);
hold on;
scatter(XYZMeasured(:,1),XYZMeasured(:,3));
scatter(XYZ(:,1),XYZ(:,3));
drawnow;
hold off; 
subplot(2,2,3);
hold on;
scatter(XYZMeasured(:,2),XYZMeasured(:,3));
scatter(XYZ(:,2),XYZ(:,3));
drawnow;
legend('truth','estimated');
hold off; 

%% Aligning point cloud by decimation



%% Aligning the point cloud (buggy)
T=[1 0 0 0; 0 1 0 0;0 0 1 0; median(XYZMeasured,1)-median(XYZ,1) 1];
tformAvg=affine3d(T);
XYZTransAvg=pctransform(pointCloud(XYZ),tformAvg);
XYZTransAvg=XYZTransAvg.Location;
close all;
XYZDecim=pcregdecim(XYZTransAvg,XYZMeasured,20,50,1);
XYZDecim=pcregdecim(XYZDecim,XYZMeasured,20,50,2);
XYZDecim=pcregdecim(XYZDecim,XYZMeasured,20,50,3);

XYZTrans=XYZDecim;
%%
    
% [tform,XYZTrans] = pcregrigid(pointCloud(XYZ),pointCloud(XYZMeasured), ...
%     'Verbose',true,'InlierRatio',0.5,'Extrapolate',true, 'InitialTransform',tformAvg);

%%
% [tform,XYZTrans] = pcregisterndt(pointCloud(XYZ),pointCloud(XYZMeasured),1000);
% XYZTrans=XYZTrans.Location;

% close all
% figure();
% subplot(2,2,1);
% scatter(XYZMeasured(:,1),XYZMeasured(:,2));
% hold on;
% scatter(XYZTransAvg(:,1),XYZTransAvg(:,2));
% scatter(XYZ(:,1),XYZ(:,2));
% hold off; 
% subplot(2,2,2);
% hold on;
% scatter(XYZMeasured(:,2),XYZMeasured(:,1));
% scatter(XYZTransAvg(:,2),XYZTransAvg(:,1));
% scatter(XYZ(:,2),XYZ(:,1));
% hold off; 
% subplot(2,2,3);
% hold on;
% scatter(XYZMeasured(:,1),XYZMeasured(:,3));
% scatter(XYZTransAvg(:,1),XYZTransAvg(:,3));
% scatter(XYZ(:,1),XYZ(:,3));
% hold off; 
% legend('Truth','transAvg','est')

figure();
subplot(2,2,1);
scatter(XYZMeasured(:,2),XYZMeasured(:,3));
hold on;
scatter(XYZTrans(:,2),XYZTrans(:,3));
hold off; 
subplot(2,2,2);
hold on;
scatter(XYZMeasured(:,2),XYZMeasured(:,1));
scatter(XYZTrans(:,2),XYZTrans(:,1));
hold off; 
subplot(2,2,3);
hold on;
scatter(XYZMeasured(:,1),XYZMeasured(:,3));
scatter(XYZTrans(:,1),XYZTrans(:,3));
hold off; 
legend('truth','estimatedMatched')
%%
close all;
figure();
% imshow(squeeze(max(cropTruth,[],3)),[],'Border','tight'); 
hold on;
scatter(XYZMeasured(:,1),XYZMeasured(:,2));
scatter(XYZ(:,1)-230,XYZ(:,2)-370);
drawnow;
hold off
legend('truth','estimatedMatched')

figure();
% imshow(squeeze(max(cropTruth,[],2)),[],'Border','tight'); 
hold on;
scatter(XYZMeasured(:,3),XYZMeasured(:,2));
scatter(XYZ(:,3)-82,XYZ(:,2)-370);
drawnow;
hold off
legend('truth','estimatedMatched')

%%

figure();
imshow(squeeze(max(cropTruth,[],3)),[],'Border','tight'); 
hold on;
scatter(XYZMeasured(:,1),XYZMeasured(:,2));
scatter(XYZTrans(:,1)+25,XYZTrans(:,2)-12);
drawnow;
hold off

figure();
imshow(squeeze(max(cropTruth,[],2)),[],'Border','tight'); 
hold on;
scatter(XYZMeasured(:,3),XYZMeasured(:,2));
scatter(XYZTrans(:,3)-19,XYZTrans(:,2)-12);
drawnow;
hold off

%%
figure();
imshow(squeeze(max(cropTruth,[],2)),[],'Border','tight'); 
hold on;
scatter(XYZMeasured(:,3),XYZMeasured(:,2));
scatter(XYZTrans(:,3),XYZTrans(:,2));
drawnow;
hold off

figure();
imshow(squeeze(max(cropTruth,[],3)),[],'Border','tight'); 
hold on;
scatter(XYZMeasured(:,1),XYZMeasured(:,2));
scatter(XYZTrans(:,1),XYZTrans(:,2));
drawnow;
hold off


%% Jaccard Index
close all;
[idx1, idx2] = colocalizationLAP(XYZMeasured, XYZTrans, 5);
JI=numel(idx1)/(size(XYZMeasured,1)+size(XYZTrans,1)-numel(idx1));
disp(['Jaccard Indx is : ' num2str(JI)]);
figure();
subplot(2,2,1);
mip=squeeze(max(cropTruth,[],3));
imshow(mip,prctile(mip(:),[5 99])); 
hold on;
scatter(XYZMeasured(:,1),XYZMeasured(:,2));
scatter(XYZTrans(:,1),XYZTrans(:,2));
scatter(XYZTrans(idx2,1),XYZTrans(idx2,2));
xlabel('X');
ylabel('Y');
hold off;

subplot(2,2,2)
mip=squeeze(max(cropTruth,[],2));
imshow(mip,prctile(mip(:),[5 99])); 
hold on;
scatter(XYZMeasured(:,3),XYZMeasured(:,2));
scatter(XYZTrans(:,3),XYZTrans(:,2));
scatter(XYZTrans(idx2,3),XYZTrans(idx2,2));
xlabel('Z');
ylabel('Y');
legend('Ground truth','False positives','True positives');
hold off;
subplot(2,2,3)
mip=squeeze(max(cropTruth,[],1))';
imshow(mip,prctile(mip(:),[5 99])); 
hold on;
scatter(XYZMeasured(:,1),XYZMeasured(:,3));
scatter(XYZTrans(:,1),XYZTrans(:,3));
scatter(XYZTrans(idx2,1),XYZTrans(idx2,3));
xlabel('X');
ylabel('Z');
hold off;

%% RMSE 
figure
histogram(sum(((XYZMeasured(idx1,:)-XYZTrans(idx2,:)).^2),2).^.5);
xlabel('distance')
RMSE=sqrt(mean(  sum(((XYZMeasured(idx1,:)-XYZTrans(idx2,:)).^2),2) ));
disp(['RMSE is : ' num2str(RMSE)]);
