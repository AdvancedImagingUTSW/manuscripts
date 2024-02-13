%% This code is to calculate the curve using max projection image of the side view
clc;
clear;
close all

acq_number=3;

imagePath = 'F:\Bingying\Curvedprops_mode\embyro1'; % load the side-view image
imageName = 'MAX_Reslice of stack1.tif';
curvename=['scanning_curve_0901-' num2str(acq_number) '.mat']; % saving file name


% the scan range and step size are
scanrange = [-150,150]; % unit um
scan_step = 0.5;
pixel_size = 150/1000; % unit um

filepath=fullfile(imagePath,imageName);


% parameters to adjust the curve
if_adjust_curve_ratio = 0;

% read the image

            InfoImage=imfinfo(filepath);
            mImage=InfoImage(1).Height;
            nImage=InfoImage(1).Width;
            NumberImages=length(InfoImage);
            
            FinalImage=zeros(mImage,nImage,NumberImages,'uint16');
            
            TifLink = Tiff(filepath, 'r');
            for i=1:NumberImages
                TifLink.setDirectory(i);
                FinalImage(:,:,i)=TifLink.read();
            end
            TifLink.close();
            
            
% segmentation 

% 1.resize the image to get isotropic sampling, and shear the image. 
resized_nImage=round((nImage-1)*scan_step/pixel_size);
resized_scanstep= pixel_size;

I = imresize(FinalImage,[mImage resized_nImage]);

tform = maketform("affine",[1 0 0; -cosd(45) sind(45) 0; 0 0 1]);
I_  = imtransform(I,tform);

% 2. perform multilevel otsu thresholding. using multithresh
level = multithresh(I_, 4); % use 4 levels 
level_seg = imquantize(I_,level);

% rebinarize
level_seg_binary = level_seg == 4;

% largest component
CC = bwconncomp(level_seg_binary);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);

level_seg_binary_largest = zeros(size(level_seg_binary,1), size(level_seg_binary,2));
level_seg_binary_largest(CC.PixelIdxList{idx})=1;

% extract midline;
level_seg_binary_largest = level_seg_binary_largest>0;
se = strel('disk',10); % size of hole to fill
level_seg_binary_largest = imclose(level_seg_binary_largest,se);
level_seg_binary_largest = imgaussfilt(level_seg_binary_largest*1., 35); % smooth the data
level_seg_binary_largest = level_seg_binary_largest/max(level_seg_binary_largest(:));
level_seg_binary_largest = level_seg_binary_largest > 0.1;
midline = bwskel(level_seg_binary_largest);

% 4. extract contour
[B,L] = bwboundaries(midline,'noholes');

scancurve = B{1};
plot(scancurve(:,1))


scancurve = scancurve(1:floor(length(scancurve)/2),:);

%% check if the ratio need to be adjusted 
if if_adjust_curve_ratio
    slope_ratio_left = 2; %percentage
    slope_ratio_right = 3; %percentage
    shift_downwards = 10; %pixels
    
    X=scancurve(:,2);
    Y=scancurve(:,1);
    [M,peak] = max(Y)
    X(peak:end)=(X(peak:end)-X(peak))*(slope_ratio_right/100+1)+X(peak);
    X(1:peak-1)=(X(1:peak-1)-X(peak-1))*(slope_ratio_left/100+1)+X(peak-1);
    
    Y = Y + shift_downwards;

    figure(100)
    imagesc(I_)
    hold on 
    plot(X, Y, 'w', 'LineWidth', 1);
    hold off 
end

%% 5.interpolate the curve based on the curve length, evenly space

if if_adjust_curve_ratio
    scancurve(:,2)= X;
    scancurve(:,1)= Y;
end

dLen=[0;vecnorm(diff( [scancurve(:,2),scancurve(:,1)] ),2,2)];

for i = 1:length(scancurve)
    Len(i) = sum(dLen(1:i));
end

l=linspace(Len(1),Len(end),round(Len(end)));
XX=interp1(Len,scancurve(:,2),l,'spline');  % use whatever method suits you the best
YY=interp1(Len,scancurve(:,1),l,'spline');  % use whatever method suits you the best

% get the ratio to the scanrange, [-150,150]shoud coorespond to size(I,2)
scan_ratio = Len(end)/size(I,2);

%%
figure(1)
imagesc(I_)
hold on 
plot(XX, YY, 'k', 'LineWidth', 3);
hold off 
saveas(gcf, fullfile(imagePath,['YZ_deskewed-' num2str(acq_number) '.png']));

%
shearing_galvo=YY/sind(45);
scanning_galvo= XX-(size(I_,1)-YY);
figure(2)
imagesc(I)
hold on
plot(scanning_galvo,shearing_galvo, 'k', 'LineWidth', 3);
saveas(gcf, fullfile(imagePath,['YZ-' num2str(acq_number) '.png']));


%%
shearing_curve=shearing_galvo-mImage/2;
% match scanning_galvo with scanning range
scanning_curve=scanning_galvo*resized_scanstep+scanrange(1);

% save_files
filepath_curve=fullfile(imagePath,curvename);
save(filepath_curve,'scanning_curve', 'shearing_curve','pixel_size','scanrange','scan_ratio');





