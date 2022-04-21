function twoStepDeconvolution(imageDirectory,imageName,psfPath,numberIterations)
%% twoStepDeconvolution
% imagePath is the path to the folder containing the image file.
% imageName is the name of the image to deconvolve
% psfPath is the path to the PSF file - e.g., /project/cellbiology/Dean_lab/shared/psfs/ctASLM2-510nm.tif
% numberIterations is usually set to 10.
%
% Written by Bo-Jui Chang, 2019.  Verified on Matlab/2019a.
%%
outputDirectory=fullfile(imageDirectory,strcat('decon_',num2str(numberIterations))); 
mkdirRobust(outputDirectory);
disp(['Data Exporting to ' outputDirectory]);

% load PSF
PSF = double(tiffRead(psfPath));

% Threshold PSF by bottom 5% & Normalize
intensityDistribution = sort(PSF(:));
PSFbackground = mean(intensityDistribution(1:size(PSF(:))/20));
disp(['The Background Intensity for the PSF is ' num2str(PSFbackground)]);
%PSF=abs(PSF-PSFbackground);
PSF=abs(PSF-118);

% Deconvolve the PSF to get a better estimate of the real PSF.
% Load the data.
filepath=fullfile(imageDirectory,imageName);
imageInfo = imfinfo(filepath);
imData = tiffRead(filepath); disp([imageName ' Loaded']);
paddedImData=padarray(single(imData),[20 20 20],'symmetric'); 
imDataMaxIntensity=max(paddedImData(:));
disp('Deconvolving Data');
[~,enhancedPSF]=deconvblind(paddedImData,PSF,numberIterations);

% Save the PSF
disp('Saving the PSF');
enhancedPSF=enhancedPSF./max(enhancedPSF(:));
enhancedPSF=uint16(enhancedPSF.*2^16);
tiffWrite(enhancedPSF,fullfile(outputDirectory,'enhancedPSF.tif'));

%% Deconvolve the Data With The Improved PSF Estimate
disp('Deconvolving Data with Enhanced PSF');
[deconvolvedImData,~]=deconvblind(paddedImData,enhancedPSF,numberIterations);
deconvolvedImData=deconvolvedImData(21:20+imageInfo(1).Height,21:20+imageInfo(1).Width,21:20+length(imageInfo));
deconvolvedImData=deconvolvedImData./max(deconvolvedImData(:));
deconvolvedImData=uint16(deconvolvedImData*imDataMaxIntensity);

% save the deconvolved image
tiffWrite(deconvolvedImData,fullfile(outputDirectory,imageName));
disp([imageName ' Deconvolution Complete']);

