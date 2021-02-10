clc; clear all;
tic
%% Input section
imagePath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/Lattice/Data_Gaussian/NA0.55,0.52/200322';
file_interest = 'Cell'; %'Squ_0,93,-93_1.00_'
number_interest = [5:12];  %%[1,2,3-5]
ChannelstoProcess= [1];
%y_ROI {1} = [1151:1200];
%y_ROI {2} = [1201:1250];
%y_ROI {3} = [1251:1300];
%y_ROI {4} = [1301:1350];
%y_ROI {5} = [1351:1400];
%y_ROI {6} = [1401:1450];
y_ROI {1} = [151:200];
y_ROI {2} = [201:250];
y_ROI {3} = [251:300];
y_ROI {4} = [301:350];
y_ROI {5} = [351:400];
y_ROI {6} = [401:450];

%pixelsToAverage_xy = 100;
%pixelsToAverage_xz = size(y_ROI{1},2);
fitRange=100; % Laterally, pixels above and below the center laterally,
fitSize=100; % Laterally, pixels to fit around the main peak, above and below the center, 3 for Bessel, 7 for HexLattice, 100 for SquLattice, 4 for 0 order
fitEdge=10; % Laterally, in pixels,
%pixelsToAverage_xz = 10;
fitRangeAxial=61; % Axially, pixels above and below the center.
fitRayleigh=45; % Axially, pixels above and below the center for Rayleigh range, need to select carefully depending on the data
cut=2; %numbers of pixels to cut in axial profile because when the light sheet is very tilted, we need to cut some black pixels. Add line 411 on 201801031
xyPixel= 81.25; %% nm
zPixel = 812.5; %% nm
%fitopt=struct('Lower',[50,1,1,1,-5,1,1], ...
%         'StartPoint',[100,300,257,350,0,10,10], ...
%              'Upper',[3000,20000,512,700,5,20,200]); %% offset, amplitude, centroid X, centroid Y, angle, width X, width Y

%% simple calculation
factor=zPixel/xyPixel;
dataNumber=1;

%% create folders
dir_SUMs=fullfile(imagePath,'SUMs');
mkdir(dir_SUMs);
mkdir(fullfile(dir_SUMs,'XY'))
mkdir(fullfile(dir_SUMs,'XY_focus'))
mkdir(fullfile(dir_SUMs,'XZ'))
mkdir(fullfile(dir_SUMs,'YZ'))
mkdir(fullfile(dir_SUMs,'analysis'))

%% MIPs

if size(number_interest,2)==0
    n = 1;
else
    n = size(number_interest,2);
end

result=['name'];
for r=1:size(y_ROI,2)
    title1=strcat('thickness_',num2str(min(y_ROI{r})),'-',num2str(max(y_ROI{r})));
    title2=strcat('length_',num2str(min(y_ROI{r})),'-',num2str(max(y_ROI{r})));
    title3=strcat('confocal parameter_',num2str(min(y_ROI{r})),'-',num2str(max(y_ROI{r})));
    result=horzcat(result,{title1,title2,title3});
end

for n1=1:n
    if size(number_interest,2)==0
        names1 = file_interest;
    else
    names1 = strcat(file_interest,num2str(number_interest(n1)));
    end
        
    
    mkdir(fullfile(dir_SUMs,'XY',names1));
    mkdir(fullfile(dir_SUMs,'XZ',names1));
    mkdir(fullfile(dir_SUMs,'YZ',names1));
    mkdir(fullfile(dir_SUMs,'three',names1));
    mkdir(fullfile(dir_SUMs,'XY_focus',names1));
    %mkdir(fullfile(dir_SUMs,'analysis',names1));
    analysisname = strcat(names1,'_fitSize',num2str(fitSize),'_fitRayleigh',num2str(fitRayleigh));
    mkdir(fullfile(dir_SUMs,'analysis',analysisname));
    %names2 = dir(fullfile(imagePath,names1));
       
    %% save parameters
fileID=fopen([fullfile(dir_SUMs,'analysis',analysisname),'/note.txt'],'wt');
fprintf(fileID,'%s%s%s\r\n','xyPixel= ',num2str(xyPixel),' nm');
fprintf(fileID,'%s%s%s\r\n','zPixelF= ',num2str(zPixel),' nm');
fprintf(fileID,'%s%s\r\n','fitRange= ',num2str(fitRange));
fprintf(fileID,'%s%s\r\n','fitSize= ',num2str(fitSize));
fprintf(fileID,'%s%s\r\n','fitEdge= ',num2str(fitEdge));
fprintf(fileID,'%s%s\r\n','fitRangeAxial= ',num2str(fitRangeAxial));
fprintf(fileID,'%s%s\r\n','fitRayleigh= ',num2str(fitRayleigh));
fprintf(fileID,'%s%s\r\n','cut= ',num2str(cut));

    filename=strcat('1_CH',num2str(ChannelstoProcess,'%02.0f'),'_000000.tif');
   
       tic
       filepath=fullfile(imagePath,names1,filename);
       %filepath=fullfile(imagePath,names1,names2(3).name);
       %imageFile = double(zeros(imageInfo(1).Width,imageInfo(1).Height,length(imageInfo)));
       InfoImage=imfinfo(filepath);
       mImage=InfoImage(1).Height;
       nImage=InfoImage(1).Width;
       NumberImages=length(InfoImage);
       
       FinalImage=zeros(mImage,nImage,NumberImages,'uint16');
       [vx,vy]=meshgrid(1:1/factor:NumberImages,1:nImage);
       [vx2,vy2]=meshgrid(1:1/factor:NumberImages,1:mImage);
       three=zeros(mImage+size(vx,2),nImage+size(vy,2));
       clear SUM*
               
    TifLink = Tiff(filepath, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   FinalImage(:,:,i)=TifLink.read();
end
TifLink.close();     

FinalImage=double(FinalImage);

[M,I] = max(max(max(FinalImage)));

center=FinalImage(:,:,I);
imwrite(uint16(center),fullfile(dir_SUMs,'XY_focus',names1,strcat('XY_',num2str(I),'_',filename)));
%imwrite(uint16(center),fullfile(dir_SUMs,'XY_focus',names1,strcat('XY_',num2str(I),'_',names2(3).name)));


sumxy=zeros(mImage,nImage);
for i=1:NumberImages
    sumxy(:,:)=sumxy(:,:)+FinalImage(:,:,i);
end
sumxy=sumxy./NumberImages;

sumxz=zeros(mImage,NumberImages);
for i=1:mImage
    tempxz(:,:)=FinalImage(i,:,:);
    sumxz(:,:)=sumxz(:,:)+tempxz(:,:);
end
sumxz=sumxz./mImage;

sumyz=zeros(nImage,NumberImages);
for i=1:nImage
    tempyz(:,:)=FinalImage(:,i,:);
    sumyz(:,:)=sumyz+tempyz(:,:);
end
sumyz=sumyz./nImage;

sumxz2=flipud(rot90(interp2(sumxz,vx,vy)));
sumyz2=interp2(sumyz,vx2,vy2);

three(1:mImage,1:nImage)=sumxy;
three(mImage+1:size(three,1),1:nImage)=sumxz2;
three(1:mImage,nImage+1:size(three,2))=sumyz2;

imwrite(uint16(sumxy),fullfile(dir_SUMs,'XY',names1,strcat('XY_',filename)));
imwrite(uint16(sumxz2),fullfile(dir_SUMs,'XZ',names1,strcat('XZ_',filename)));
imwrite(uint16(sumyz2),fullfile(dir_SUMs,'YZ',names1,strcat('YZ_',filename)));
imwrite(uint16(three),fullfile(dir_SUMs,'three',names1,strcat('three_',filename)));
%imwrite(uint16(sumxy),fullfile(dir_SUMs,'XY',names1,strcat('XY_',names2(3).name)));
%imwrite(uint16(sumxz2),fullfile(dir_SUMs,'XZ',names1,strcat('XZ_',names2(3).name)));
%imwrite(uint16(sumyz2),fullfile(dir_SUMs,'YZ',names1,strcat('YZ_',names2(3).name)));
%imwrite(uint16(three),fullfile(dir_SUMs,'three',names1,strcat('three_',names2(3).name)));

result_temp=[];
lateral_raw=[];
axial_raw=[];
confocal_raw=[];

%% Preparing for rotating the stack in xy (laterally)
   %% Identify in-focus image-plane in Z.
        planeMaxs = zeros(size(FinalImage,3),1);
        planeNumber = 1:1:size(FinalImage,3);
        parfor planeIdx = 5:1:size(FinalImage,3)-5
            planeMaxs(planeIdx) = max(max(FinalImage(:,:,planeIdx)));
        end
        
        % Fit model to data.
        [xData, yData] = prepareCurveData( planeNumber, planeMaxs );
        ft = fittype( 'a1*exp(-((x-b1)/c1)^2)+d');
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off'; opts.Lower = [0 0 1 0];
        opts.StartPoint = [max(planeMaxs) size(FinalImage,3)/2 10 planeMaxs(1)];
        opts.Upper = [max(planeMaxs)*2 Inf Inf max(planeMaxs)];
        [fitresult, gof] = fit( xData, yData, ft, opts );
        inFocusIdx = round(fitresult.b1);

    %% Analyze In-Focus Image
        % Isolate the Beam Waist
        %inFocusPlane = imageVol(:,:,inFocusIdx);
        inFocusPlane = FinalImage(:,:,inFocusIdx);
        figure(1); imagesc(inFocusPlane);
        % Background Subtract
        inFocusPlane = inFocusPlane-min(inFocusPlane(:));
        
        parfor i = 1:1:size(inFocusPlane,1);
            [locs(i), pks(i)] = max(inFocusPlane(i,:));
        end
        
        % Fit model to data.
        [xData, yData] = prepareCurveData([1:size(inFocusPlane,1)], double(pks));
        ft = fittype( 'a1*x+b1', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [-Inf 0];
        opts.StartPoint = [0 mean(pks)];
        opts.Upper = [Inf Inf];
        [fitresult, gof] = fit( xData, yData, ft, opts );
        inFocusIdx = round(fitresult.b1);
        rotAngleXY=atand(-fitresult.a1);
        % Rotate the Data To Account for Non-Perfect Alignment
        inFocusRotatePlane = imrotate(inFocusPlane,rotAngleXY,'bicubic');
        figure(1); imagesc(inFocusRotatePlane);
        disp(['Real rotAngleXY = ' num2str(rotAngleXY), ' degree']);
%         subplot(1,3,1);
%         h = plot( fitresult, xData, yData );
%         legend( h, 'Beam Rotation', 'Linear Fit', 'Location', 'NorthEast' );
%         subplot(1,3,2:3); 
%         imshowpair(inFocusPlane,inFocusRotatePlane,'montage');

%% Preparing for rotating the stack in xz direction (axially) 
    % use the whole stack to rotate a little in xz direction.
    % 1. average the profile in xz direction.
    % 2. select a range (fitRangeAxial). 
    % 3. fit Gaussian curve to each z position
    % 4. use the center of each fitted curve to determine the rotation angle.
        close all
        FinalImage2=squeeze(mean(FinalImage,1));
        [pks1, locs1] = max(FinalImage2); 
        [M,I]=max(pks1);
        %FinalImage2=FinalImage2(:,(I-fitRangeAxial):(I+fitRangeAxial));
        center=floor(size(FinalImage2,2)/2)+1;
        FinalImage2=FinalImage2(:,(center-fitRangeAxial):(center+fitRangeAxial));
        figure(2);imagesc(FinalImage2);
      
        xData=(1:size(FinalImage2,1));
        
        for planeIdx = 1:size(FinalImage2,2)
            [xData, yData] = prepareCurveData(xData,FinalImage2(:,planeIdx));
            ft = fittype( 'a1*exp(-((x-b1)/c1)^2)+d', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            %opts.Lower = [10 0 0 0];
            opts.StartPoint = [max(FinalImage2(:,planeIdx)) size(FinalImage2,1)/2 20 min(FinalImage2(:,planeIdx))];
            %opts.Upper = [max(planeMaxs)*2 Inf Inf max(planeMaxs)];
            [fitresult, gof] = fit( xData, yData, ft, opts );
            b1(planeIdx)=fitresult.b1;
        end
                 
        [xData, yData] = prepareCurveData(1:size(b1,2), b1(:));
        ft = fittype( 'a1*x+b1', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [-Inf 0];
        opts.StartPoint = [0 mean(b1(:))];
        opts.Upper = [Inf Inf];
        [fitresult, gof] = fit( xData, yData, ft, opts );
        figure(3);plot(fitresult,xData,yData);
        rotAngleXZ=atand(fitresult.a1);
        FinalImage2Rotate = imrotate(FinalImage2,rotAngleXZ,'bicubic','crop');
        figure(2); imagesc(FinalImage2Rotate);
        disp(['Real rotAngleXZ = ' num2str(atand(fitresult.a1/factor)), ' degree']);
%% Measure the light-sheet properties in each ROI

for r=1:size(y_ROI,2)
sumxz_ROI=zeros(mImage,NumberImages);
xz_ROI=y_ROI{r};
pixelsToAverage_xz = size(y_ROI{r},2);
for i=min(xz_ROI):max(xz_ROI)
    tempxz_ROI(:,:)=FinalImage(i,:,:);
    sumxz_ROI(:,:)=sumxz_ROI(:,:)+tempxz_ROI(:,:);
end
sumxz_ROI=sumxz_ROI./size(xz_ROI,2);
sumxz_ROI2=flipud(rot90(interp2(sumxz_ROI,vx,vy)));
    
newname=strcat('_',num2str(min(xz_ROI)),'-',num2str(max(xz_ROI)));
mkdir(fullfile(dir_SUMs,strcat('XZ',newname),names1));
imwrite(uint16(sumxz_ROI2),fullfile(dir_SUMs,strcat('XZ',newname),names1,strcat('XZ_',filename)));
%imwrite(uint16(sumxz_ROI2),fullfile(dir_SUMs,strcat('XZ',newname),names1,strcat('XZ_',names2(3).name)));

        
    %% Laterally Average beam to Remove Lattice Structure.
        close all
        
        inFocusPlane = inFocusRotatePlane;
        center_ROI= xz_ROI(ceil(end/2));
        
        %pixelsToAverage=10;
        %beamLateralCrossSection = inFocusPlane(size(inFocusPlane,2)/2-round(pixelsToAverage_xy/2):size(inFocusPlane,2)/2+round(pixelsToAverage_xy/2),:);
        %beamLateralCrossSection = squeeze(sum(beamLateralCrossSection,1));
        %beamLateralCrossSection = inFocusPlane(center_ROI-round(pixelsToAverage_xy/2):center_ROI+round(pixelsToAverage_xy/2),:);
        beamLateralCrossSection = inFocusPlane(center_ROI-round(pixelsToAverage_xz/2):center_ROI+round(pixelsToAverage_xz/2),:);
        beamLateralCrossSection = squeeze(mean(beamLateralCrossSection,1)); %change from sum to mean @20200401
        
        % Truncate the Data to Remove Large Region Beyond Beam Focus
        [~, locs] = max(beamLateralCrossSection);
        %fitSize = 100;
        %fitSize = pixelsToAverage_xz;
        beamLateralCrossSection = beamLateralCrossSection(locs-fitRange:locs+fitRange);
        %beamLateralCrossSection = beamLateralCrossSection(locs-fitSize:locs+fitSize);
        
        % Prepare Curve Fit Data, Fittype, and Options
        pixelNumber = 1:1:size(beamLateralCrossSection,2);
        xAxis = pixelNumber*xyPixel/1000;
        
        [xData, yData] = prepareCurveData( xAxis, beamLateralCrossSection );
        ft = fittype( 'a1*exp(-((x-b1)/c1)^2)+d');
        %ft = fittype( 'a1*exp(-((x-b1)/c1)^2)+d', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        startps = [max(beamLateralCrossSection) xData(ceil(end/2),1) (fitSize/2)*xyPixel/1000 mean(beamLateralCrossSection(1:fitEdge))];
        exclx= ((xData > fitEdge*xyPixel/1000) & (xData < (ceil(size(pixelNumber,2)/2)-fitSize)*xyPixel/1000)) | ((xData > (ceil(size(pixelNumber,2)/2)+fitSize)*xyPixel/1000) & (xData < (length(beamLateralCrossSection)-fitEdge)*xyPixel/1000)) ;
        %opts.Lower = [0 0 0 0];
        %opts.StartPoint = [max(beamLateralCrossSection) xData(ceil(end/2),1) (fitSize/2)*xyPixel/1000 mean(beamLateralCrossSection(1:fitEdge))]
        %opts.Upper = [max(beamLateralCrossSection)*2 Inf Inf max(beamLateralCrossSection)];
        %opts.Exclude= ((xData > fitEdge*xyPixel/1000) & (xData < (locs-fitSize)*xyPixel/1000)) | ((xData > (locs+fitSize)*xyPixel/1000) & (xData < (length(beamLateralCrossSection)-fitEdge)*xyPixel/1000)) ;
        % fit(xData, yData,ft,)
        % Fit model to data.
        fitresult = fit( xData, yData, ft, 'Start', startps, 'Exclude', exclx );
        %[fitresult, gof] = fit( xData, yData, ft, opts );
        %[fitresult, gof] = fit( xData, yData, ftlat, 'Start', StartPoint, 'Exclude', exclx );
        % Calculate Beam Size
        beamCenter=fitresult.b1;
        beamLateralFWHM = 2*sqrt(log(2))*fitresult.c1;
        disp(['Beam Lateral FWHM in Microns = ' num2str(beamLateralFWHM)])
        
        % Beam Confinment
        % 15 pixels is 975 nm
        [pks, locs] = max(beamLateralCrossSection);
        inFocus = sum(beamLateralCrossSection(ceil(size(pixelNumber,2)/2)-fitSize:ceil(size(pixelNumber,2)/2)+fitSize));
        %inFocus = sum(beamLateralCrossSection(locs-fitSize:locs+fitSize));
        totalIntensity = sum(beamLateralCrossSection);
        percentInFocus = (inFocus./totalIntensity)*100;
        
        % Plot fit with data.
        figure(2)
        subplot(1,2,1); 
        imshow(inFocusPlane,[]);
        
        subplot(1,2,2); 
        h = plot( fitresult, xData, yData );
        legend( h, 'Beam Cross-Section', 'Gaussian Fit', 'Location', 'NorthWest' );
        xlabel Microns; 
        ylabel Intensity; 
        grid on; 
        xlim([0 max(xData)]);
        text(5,max(yData),['Beam FWHM = ' num2str(beamLateralFWHM) '\mum'])
        text(5,max(yData)-max(yData)/10,['% Confinement = ' num2str(percentInFocus)]);
        legend('off');
        %print(['Cell' num2str(dataNumber) '_lateral'],'-depsc','-tiff')
        %savefig(['Cell' num2str(dataNumber) '_lateral']);
        
        figure(3)
        h = plot( fitresult, xData, yData );
        legend( h, 'Beam Cross-Section', 'Gaussian Fit', 'Location', 'NorthWest' );
        xlabel Microns; 
        ylabel Intensity; 
        grid on; 
        xlim([0 max(xData)]);
        text(0,max(yData),['Beam FWHM = ' num2str(beamLateralFWHM) '\mum'])
        text(0,max(yData)-max(yData)/10,['% Confinement = ' num2str(percentInFocus)]);
        legend('off');
        title(['Beam FWHM = ' num2str(beamLateralFWHM),' um']);
                
        %saveas(figure(3),fullfile(dir_SUMs,'analysis',names1,strcat('xy_',names2(3).name)));
        %saveas(figure(3),fullfile(dir_SUMs,'analysis',names1,strcat('xy',newname,'.tif')));
        saveas(figure(3),fullfile(dir_SUMs,'analysis',analysisname,strcat('xy',newname,'.tif')));
        
        lateral_raw(1:length(xData),2*r-1)=xData;  
        lateral_raw(1:length(yData),2*r)=yData; 
         %%
        % for Cell 6, should be exactly 10 pixels apart.  10Hz?
        %close all
        
        dataFFT = fft(yData,size(yData,1));
        Pyy = dataFFT.*conj(dataFFT)/size(yData,1);
        f = 1/size(yData,1)*(0:floor(size(yData,1)./2));
        figure (1)
        plot(f,Pyy(1:round(size(yData,1)./2)))
        %f = 103, Pyy = 101.
        title('Power spectral density')
        xlabel('Frequency (Hz)')
        set(gca, 'XScale', 'log')
        %print(['Cell' num2str(dataNumber) '_fft'],'-depsc','-tiff')
        %savefig(['Cell' num2str(dataNumber) '_fft']);
                   
        %% Measure Light-sheet axial properties
        % Want to measure how thick the beam is in Z, a the brighest
        % position of the in-focus image.
        
        subFinalImage=FinalImage(xz_ROI,:,:);
        subFinalImage=squeeze(mean(subFinalImage,1));
        subFinalImage=subFinalImage(:,(center-fitRangeAxial):(center+fitRangeAxial));
        %subFinalImage=subFinalImage(:,(I-fitRangeAxial):(I+fitRangeAxial));
        figure(1);
        subplot(2,1,1);imagesc(subFinalImage);
        subFinalImage = imrotate(subFinalImage,rotAngleXZ,'bicubic','crop');   
        subplot(2,1,2);imagesc(subFinalImage);
        
        %% Measure Rayleigh length (confocal parameter)        
        
        [pks1, locs1] = max(subFinalImage); 
        [M,I]=max(pks1);
        if I-(fitRayleigh+10) <= 1 | I+(fitRayleigh+10) >= size(subFinalImage,2)
            lowerend = 1 
            upperend = size(subFinalImage,2)
            %lowerend = round(length(pks1)/2)-(fitRayleigh+10)
            %upperend = round(length(pks1)/2)+(fitRayleigh+10)
        else
            lowerend = I-(fitRayleigh+10)
            upperend = I+(fitRayleigh+10)
        end
        
        %if I-(fitRayleigh+10) <= 1
        %    lowerend = 1
        %else
        %    lowerend = I-(fitRayleigh+10)
        %end
        %if I+(fitRayleigh+10) >= size(subFinalImage,2)
        %    upperend = size(subFinalImage,2)
        %else
        %    upperend = I+(fitRayleigh+10)
        %end
            
        CPImage=subFinalImage(round(mean(locs1))-(fitSize+5):round(mean(locs1))+(fitSize+5),lowerend:upperend);
        %CPImage=subFinalImage(:,(I-(fitRayleigh+5)):(I+(fitRayleigh+5)));
         
        xData= 1:size(CPImage,1);
               

        for planeIdx = 1:size(CPImage,2)
        [xData, yData] = prepareCurveData(xData,CPImage(:,planeIdx));
        ft = fittype( 'a1*exp(-((x-b1)/c1)^2)+d', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            %opts.Lower = [10 0 0 0];
            opts.StartPoint = [max(CPImage(:,planeIdx)) size(CPImage,1)/2 20 min(CPImage(:))];
            %opts.Upper = [max(planeMaxs)*2 Inf Inf max(planeMaxs)];
            [fitresult, gof] = fit( xData, yData, ft, opts );
            rayleighSigma(planeIdx) = (xyPixel*2*sqrt(log(2))*fitresult.c1)/1000;
            %rayleighSigma(planeIdx) = (xyPixel*sqrt(2*log(2))*fitresult.c1)/1000;
        end
        
        rayleighX = 1:1:length(rayleighSigma);
        rayleighX = rayleighX*zPixel/1000;
        
        [M,I]=min(rayleighSigma);
        if I-fitRayleigh <= 1 | I+fitRayleigh >= length(rayleighX)
            lowerend = round(length(rayleighX)/2)-fitRayleigh
            upperend = round(length(rayleighX)/2)+fitRayleigh
        else
            lowerend = I-fitRayleigh
            upperend = I+fitRayleigh
        end
       
        %x1=rayleighX;
        x1=rayleighX(lowerend:upperend);
        %x1=rayleighX(I-fitRayleigh:I+fitRayleigh);
        %y1=rayleighSigma;
        y1=rayleighSigma(lowerend:upperend);
        %y1=rayleighSigma(I-fitRayleigh:I+fitRayleigh);
        xscale=min(x1):xyPixel/1000:max(x1);
        rescale=interp1(x1,y1,xscale,'spline');
        [xData, yData] = prepareCurveData(xscale, rescale);
        
       % Rayleigh Equation
       ft = fittype( 'a1*sqrt(1+((x-b1)/c1)^2)', 'independent', 'x', 'dependent', 'y' );
       opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
       opts.Display = 'Off';
       opts.StartPoint = [min(yData) mean(xData) 2];
       opts.Lower = [min(yData)-0.2 mean(xData)-1 0];
       opts.Upper = [min(yData)+0.2 mean(xData)+1 20];
       [fitresult, gof] = fit( xData, yData, ft, opts );
       Confocal_par= 2*fitresult.c1;
       figure(4)
       plot(fitresult,rayleighX,rayleighSigma);
       xlim([0 max(rayleighX)]);
       disp(['raw 2w0 in Microns = ' num2str(min(yData))])
       disp(['2w0 in Microns = ' num2str(fitresult.a1)])
       disp(['2wr in Microns = ' num2str(fitresult.a1*sqrt(2))])
       fprintf(fileID,'%s%s%s\r\n','raw 2w0= ',num2str(min(yData)),' um');
       fprintf(fileID,'%s%s%s\r\n','2w0= ',num2str(fitresult.a1),' um');
       fprintf(fileID,'%s%s%s\r\n','2wr= ',num2str(fitresult.a1*sqrt(2)),' um');
   %   % Gaussian curve fitting method
   %     ft = fittype( 'a1*exp(((x-b1)/c1)^2)+d', 'independent', 'x', 'dependent', 'y' );
   %     opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); 
   %     opts.Display = 'Off';
   %     opts.Lower = [-Inf mean(xData)-2 0 -Inf];
   %     opts.StartPoint = [mean(yData) mean(xData) 2 min(yData)];
   %     opts.Upper = [Inf mean(xData)+2 4 Inf];
   %     [fitresult, gof] = fit( xData, yData, ft, opts );
   %  
   %     Rlength = 2*sqrt(log(2))*fitresult.c1; 
   %     figure (4)
   %     plot( fitresult, xData, yData )
        
   %    % parabolic curve fitting method       
   %     p=polyfit(x1,y1,2); %fit 2nd degree polynomials
   %     %y2=polyval(p,rayleighX);
   %     equ=poly2sym(p);
   %             
   %     fun=matlabFunction(equ);
   %     [s0, p0] = fminbnd(fun,min(rayleighX),max(rayleighX)); %find the minimum of the fitting polynomials
   %     S=vpasolve(equ==sqrt(2)*p0);
   %     Rlength=round(single(abs(S(1)-S(2))),2); %find the Rayleigh length;
         disp(['Confocal_par. in Microns = ' num2str(Confocal_par)])
        
   %     figure (4)
   %     plot(rayleighX, rayleighSigma);
   %     hold on
   %     %plot(rayleighX,y2);
   %     %hold on
   %     fplot(equ,[min(rayleighX),max(rayleighX)]);
   %     hold off
        
        xlabel('Microns'); ylabel('Beam Waist (microns)'); grid on;
        legend('off'); 
        title(['Confocal par. = ' num2str(Confocal_par),' um']);
        %print(['Cell' num2str(dataNumber) '_rayleigh'],'-depsc','-tiff')
        %savefig(['Cell' num2str(dataNumber) '_rayleigh']);
        saveas(figure(4),fullfile(dir_SUMs,'analysis',analysisname,strcat('Confocal_par',newname,'.tif')));
        
        Rayleigh_raw(1:length(rayleighX),2*r-1)=rayleighX;  
        Rayleigh_raw(1:length(rayleighSigma),2*r)=rayleighSigma;  
        
         %% Measure Rayleigh Length (original Rayleigh length measurement, but not what I want ...)
        % Want to measure how thick the beam is in Z, a the brighest
        % position of the in-focus image.
%        close all
%        
%        inFocusImage = FinalImage(:,:,inFocusIdx);
%        [pks1, locs1] = max(inFocusImage);
%        [pks2, locs2] = max(pks1);  
%        
%        for planeIdx = 1:size(FinalImage,3)
%            [xData, yData] = prepareCurveData(1:size(FinalImage,2), double(FinalImage(locs2,:,planeIdx)));
%            ft = fittype( 'a1*exp(-((x-b1)/c1)^2)+d', 'independent', 'x', 'dependent', 'y' );
%            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%            opts.Display = 'Off';
%            opts.Lower = [10 0 0 0];
%            opts.StartPoint = [max(max(FinalImage(:,:,planeIdx))) size(FinalImage,2)/2 40 min(FinalImage(:))];
%            opts.Upper = [max(planeMaxs)*2 Inf Inf max(planeMaxs)];
%            [fitresult, gof] = fit( xData, yData, ft, opts );
%            rayleighSigma(planeIdx) = (xyPixel*2*sqrt(log(2))*fitresult.c1)/1000;
%            %rayleighSigma(planeIdx) = (xyPixel*sqrt(2*log(2))*fitresult.c1)/1000;
%        end
%
%        rayleighX = 1:1:length(rayleighSigma);
%        rayleighX = rayleighX*zPixel/1000;
%        figure (1)
%        plot(rayleighX, rayleighSigma);
%        xlabel('Beam Length (microns)'); ylabel('Beam Waist (microns)'); grid on;
%        legend('off'); ylim([0 10]);
%        title(['Rayleigh = ' num2str(rayleighSigma(planeIdx))]);
%        %print(['Cell' num2str(dataNumber) '_rayleigh'],'-depsc','-tiff')
%        %savefig(['Cell' num2str(dataNumber) '_rayleigh']);
%        saveas(figure(1),fullfile(dir_SUMs,'analysis',strcat(names1,'_fitSize',num2str(fitSize)),strcat('Rayleigh',newname,'.tif')));
         
               
                %close all
        
        % Find Peak
        
        %% Beam Propagation Length
        
        [maxAxial, locsAxial] = max(subFinalImage');
        [~, locsAxialPeak] = max(maxAxial);
        
        % Remove Line Profile
        beamPropagation = subFinalImage';
        beamPropagation = beamPropagation(:,locsAxialPeak);
        beamPropagation = beamPropagation(cut+1:size(beamPropagation,1)-cut); % when the light sheet is very tilted, we need to cut some black pixels
        
        
        % Prepare x-Axis in Microns
        xAxis = (zPixel/1000)*[1:size(beamPropagation,1)];
        
        % Fit model to data.  
        [xData, yData] = prepareCurveData(xAxis, beamPropagation');
        ft = fittype( 'a1*exp(-((x-b1)/c1)^2)+d', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); 
        opts.Display = 'Off';
        opts.Lower = [mean(yData) 0 0 min(yData)];
        opts.StartPoint = [max(yData)-min(yData) mean(xAxis) 5 min(yData)];
        opts.Upper = [max(beamPropagation)*2 Inf Inf max(beamPropagation)];
        [fitresult, gof] = fit( xData, yData, ft, opts );

        beamAxialFWHM = 2*sqrt(log(2))*fitresult.c1; 
        disp(['Beam Axial FWHM in Microns = ' num2str(beamAxialFWHM)])
        
        figure (5)
        subplot(1,2,1); 
        %imshow(subFinalImage,[]);
        imagesc(subFinalImage);
        
        subplot(1,2,2); 
        h = plot( fitresult, xData, yData );
        xlim([0 max(xData)]);
        xlabel Microns; 
        ylabel Intensity; 
        grid on
        legend('off');
        title(['Beam Length = ' num2str(beamAxialFWHM),' um']);
        %print(['Cell' num2str(dataNumber) '_axial'],'-depsc','-tiff')
        %savefig(['Cell' num2str(dataNumber) '_axial']);     
        
        figure(6)
        h = plot( fitresult, xData, yData );
        xlim([0 max(xData)]);
        xlabel Microns; 
        ylabel Intensity; 
        grid on
        legend('off');
        title(['Beam Length = ' num2str(beamAxialFWHM),' um']);
        
        %saveas(figure(3),fullfile(dir_SUMs,'analysis',names1,strcat('xz_',names2(3).name)));
        %saveas(figure(3),fullfile(dir_SUMs,'analysis',names1,strcat('xz',newname,'.tif')));
        saveas(figure(6),fullfile(dir_SUMs,'analysis',analysisname,strcat('xz',newname,'.tif')));
        
        axial_raw(1:length(xData),2*r-1)=xData;  
        axial_raw(1:length(yData),2*r)=yData;  
        %% Beam Propagation Length (original code, but I don't think I need it anymore because I already have the image rotated and can go straight fitting the Gaussian curve @20200402)
        %close all
        % Interpolate Z.
        % % Interpolate the Z-Axis
        
        
%        imageVolInterp = interpolateZ(FinalImage, xyPixel, zPixel);
%        disp('Data Interpolated');
%        
%        % Analyze Central Slice in X-Direction
%        %centralSlice = round(size(imageVolInterp,1)/2);
%        %pixelsToAverage = 100;
%        %axialCrossSection = squeeze(sum(imageVolInterp(centralSlice-round(pixelsToAverage_xz./2):centralSlice+round(pixelsToAverage_xz./2),:,:)));
%        axialCrossSection = squeeze(sum(imageVolInterp(center_ROI-round(pixelsToAverage_xz./2)+1:center_ROI+round(pixelsToAverage_xz./2),:,:)));
%        
%        
        %% Account for Beam Angle
%        [maxAxial, locsAxial] = max(axialCrossSection);
%        [~, locationAxialFocus] = max(maxAxial);
%        
%        numberPixelsAxial = 100; %was100 was 300
%        startPixel = locationAxialFocus-round(numberPixelsAxial./2);
%        endPixel = locationAxialFocus+round(numberPixelsAxial./2);
        
%        yData = locsAxial(startPixel:endPixel);
%        xData = startPixel:endPixel;
%        % Fit model to data.
%        [xData, yData] = prepareCurveData(xData, yData);
%        ft = fittype( 'a1*x+b1', 'independent', 'x', 'dependent', 'y' );
%        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%        opts.Display = 'Off';
%        opts.Lower = [-Inf 0];
%        opts.StartPoint = [0 mean(pks)];
%        opts.Upper = [Inf Inf];
%        [fitresult, gof] = fit( xData, yData, ft, opts );
%        inFocusIdx = round(fitresult.b1);
        
%        % Rotate the Data To Account for Non-Perfect Alignment
%        axialCrossSectionRotate = imrotate(axialCrossSection,atand(fitresult.a1),'bicubic');
%        % subplot(1,3,1); 
%        % h = plot(fitresult, xData, yData );
%        % legend( h, 'Max Intensity', 'Linear Fit', 'Location', 'NorthEast' );
%        % subplot(1,3,2:3); 
%        % imshowpair(axialCrossSection,axialCrossSectionRotate,'montage');
%         
        
        %% Analyze the Rotated Data To Determine the Beam Length
%        %close all
%        
%        % Find Peak
%        [maxAxial, locsAxial] = max(axialCrossSectionRotate');
%        [~, locsAxialPeak] = max(maxAxial);
%        
%        % Remove Line Profile
%        beamPropagation = axialCrossSectionRotate';
%        beamPropagation = beamPropagation(:,locsAxialPeak);
%        beamPropagation = beamPropagation(cut+1:size(beamPropagation,1)-cut); % when the light sheet is very tilted, we need to cut some black pixels
%             
%        % Prepare x-Axis in Microns
%        xAxis = (xyPixel/1000)*[1:size(beamPropagation,1)];
%        
%        % Fit model to data.  
%        [xData, yData] = prepareCurveData(xAxis, beamPropagation');
%        ft = fittype( 'a1*exp(-((x-b1)/c1)^2)+d', 'independent', 'x', 'dependent', 'y' );
%        opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); 
%        opts.Display = 'Off';
%        opts.Lower = [mean(yData) 0 0 min(yData)];
%        opts.StartPoint = [max(yData)-min(yData) mean(xAxis) 5 min(yData)];
%        opts.Upper = [max(beamPropagation)*2 Inf Inf max(beamPropagation)];
%        [fitresult, gof] = fit( xData, yData, ft, opts );
%
%        beamAxialFWHM = 2*sqrt(log(2))*fitresult.c1; 
%        disp(['Beam Axial FWHM in Microns = ' num2str(beamAxialFWHM), 'um'])
%        
%        figure (4)
%        subplot(1,2,1); 
%        imshow(axialCrossSectionRotate,[]);
%        
%        subplot(1,2,2); 
%        h = plot( fitresult, xData, yData );
%        xlim([0 max(xData)]);
%        xlabel Microns; 
%        ylabel Intensity; 
%        grid on
%        legend('off');
%        title(['Beam Length = ' num2str(beamAxialFWHM)]);
%        print(['Cell' num2str(dataNumber) '_axial'],'-depsc','-tiff')
%        savefig(['Cell' num2str(dataNumber) '_axial']);     
%        
%        figure(5)
%        h = plot( fitresult, xData, yData );
%        xlim([0 max(xData)]);
%        xlabel Microns; 
%        ylabel Intensity; 
%        grid on
%        legend('off');
%        title(['Beam Length = ' num2str(beamAxialFWHM)]);
        
%        %saveas(figure(3),fullfile(dir_SUMs,'analysis',names1,strcat('xz_',names2(3).name)));
%        %saveas(figure(3),fullfile(dir_SUMs,'analysis',names1,strcat('xz',newname,'.tif')));
%        saveas(figure(5),fullfile(dir_SUMs,'analysis',strcat(names1,'_fitSize',num2str(fitSize)),strcat('xz',newname,'.tif')));
            
%
%        %close (figure(3))       
   
    result_temp=horzcat(result_temp,{beamLateralFWHM,beamAxialFWHM,Confocal_par});
           
toc
end
   
fprintf(fileID,'%s%s%s\r\n','Real rotAngleXY= ',num2str(rotAngleXY),' degree');
fprintf(fileID,'%s%s%s\r\n','Real rotAngleXZ= ',num2str(rotAngleXZ),' degree');
fclose(fileID);

    lateralname=fullfile(imagePath,strcat(names1,'_lateral_raw_fitSize',num2str(fitSize),'_fitRayleigh',num2str(fitRayleigh),'.csv'));
    xlswrite(lateralname,lateral_raw);

    axialname=fullfile(imagePath,strcat(names1,'_axial_raw_fitSize',num2str(fitSize),'_fitRayleigh',num2str(fitRayleigh),'.csv'));
    xlswrite(axialname,axial_raw);

    Confocalname=fullfile(imagePath,strcat(names1,'_confocal_raw_fitSize',num2str(fitSize),'_fitRayleigh',num2str(fitRayleigh),'.csv'));
    xlswrite(Confocalname,Rayleigh_raw);
    
    result=vertcat(result,[names1,result_temp]);
    
end

%resultname=fullfile(imagePath,strcat(file_interest,num2str(number_interest),'.csv'));
resultname=fullfile(imagePath,strcat(file_interest,num2str(number_interest),'_fitSize',num2str(fitSize),'_fitRayleigh',num2str(fitRayleigh),'.csv'));
fid = fopen(resultname,'wt');
if fid>0
     for k=1:size(result,1)
         fprintf(fid,repmat('%s,',[1,size(result,2)]),result{k,:});
         fprintf(fid,'\r\n');
     end
     fclose(fid);
end
 
%    %lateralname=fullfile(imagePath,strcat(file_interest,num2str(number_interest),'_lateral_raw.csv'));
%    lateralname=fullfile(imagePath,strcat(file_interest,num2str(number_interest),'_lateral_raw_fitSize',num2str(fitSize),'.csv'));
%    xlswrite(lateralname,lateral_raw);
%    %axialname=fullfile(imagePath,strcat(file_interest,num2str(number_interest),'_axial_raw.csv'));
%    axialname=fullfile(imagePath,strcat(file_interest,num2str(number_interest),'_axial_raw_fitSize',num2str(fitSize),'.csv'));
%    xlswrite(axialname,axial_raw);
%    Confocalname=fullfile(imagePath,strcat(file_interest,num2str(number_interest),'_confocal_raw_fitSize',num2str(fitSize),'.csv'));
%    xlswrite(Confocalname,Rayleigh_raw);

toc
disp('All Done');

