%% Manuscript Information
%
% Authors: Tadamoto Isogai, et al. 
% Title: Direct Arp2/3-vinculin binding is required for pseudopod extension, 
% but only on compliant substrates and in 3D
% 
% Abstract: A critical step in cell morphogenesis is the extension of actin-dense pseudopods, 
% controlled by actin-binding proteins (ABPs). The relationship between actin regulation 
% and pseudopod extension is well-understood for cells cultured on glass coverslips, 
% but much less so in compliant three-dimensional environments. Here, we knocked-out a 
% series of ABPs in osteosarcoma cells and evaluated the effect of each ABP on pseudopod 
% extension on glass surfaces (2D) and in reconstituted collagen gels (3D). Cells 
% lacking the longest Arp3 isoform, or with attenuated Arp2/3 activity,
% had the strongest reduction in pseudopod formation in 3D environments. 
% This phenotype largely stemmed from reduced activity of the recently 
% characterized hybrid Arp2/3-vinculin complex, whereas its activity 
% was dispensable for pseudopod extension in 2D. Our data suggests that 
% the synergy of concurrent formation of actin branches and nascent 
% adhesions provided by the hybrid Arp2/3-vinculin complex is essential
% to overcome the reduced availability of mechanically stable links between
% fibrous extracellular matrix and actin cytoskeleton in 3D. This highlights
% how experiments on stiff, planar substrates may conceal actin architectural 
% features that are essential for morphogenesis in 3D.

%% Analysis of 2D spreading

close all;
clc; clear;

% Define rawData location
dataDirectory = '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/raw/tisogai/H-EPI/20211005_CK666inARP3KOs/DMSO';
disp('Analyzing following folders:')
ls(dataDirectory);

% Define Location to Save
tempPath = '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/tisogai/3D-ARP3 Project/U2OS-2D/211005_CK666inARP3KO_2D_exp1/DMSO';
analysisFolder = 'stats summary';
mkdir(tempPath,analysisFolder);
savePath = [tempPath,'/',analysisFolder];
cd(savePath);

% Define experiment id
pathDate = '';

% Define pixel size:
PhysicalSize = 0.65; % um
Small = round(60);
Big = round(10000);

%% Directory contents
directoryContents = dir(dataDirectory);
for dirIdx = 3:1:length(directoryContents)
    folderPath{dirIdx-2} = directoryContents(dirIdx).name;
end

%% Watershed segmentation and feature extraction
for folderIdx = 1:size(folderPath,2)
    disp(['Processing ' folderPath{folderIdx}]);
    folderPathLong = fullfile(dataDirectory,folderPath{folderIdx});
    subFolder = dir(folderPathLong);
    [numberCells] = numberOfStrings(folderPathLong,folderPath{folderIdx});
    
    %% Load Data
    for subfolderIdx = 3:1:length(subFolder)
        disp(['Processing Cell Number ' num2str(subfolderIdx - 2)]);
        
        try
            imagePath = fullfile(folderPathLong,subFolder(subfolderIdx).name);
        catch
            imagePath = fullfile(folderPathLong,subFolder(subfolderIdx).name);
        end
        
        filePath = dir(fullfile(imagePath, '*.tif'));
        imData = imread(fullfile(imagePath,filePath.name));
        Iadj = imadjust(imData);
        
        %% Denoise the Data
        Iadj2 = double(Iadj);
        residual=(Iadj2-medfilt2(Iadj2,[3 3]));
        sigma = std(residual(:));
        
        %Identify and Remove Bright Dots
        thresh = mean(residual(:))+3*sigma;
        indx=residual>thresh;
        filteredData=Iadj;
        f=medfilt2(Iadj,[3 3]);
        filteredData(indx) = f(indx);
        
        %Identify and Remove Dark Pixels
        thresh = mean(residual(:))-3*sigma;
        clear indx
        indx=residual<thresh;
        filteredData(indx) = f(indx);
                
        % Background substraction and Mark the Foreground objects
        se = strel('disk',30);
        background = imopen(filteredData,se);
        imDataBG = filteredData - background;
        
        % Show the Filtered Data
        figure; 
        subplot(1,3,1); imshow(Iadj)
        subplot(1,3,2); imshow(filteredData,[]);
        subplot(1,3,3); imshow(imDataBG,[]);
        
        %% Threshold Data - Version 1
        imDataBG_BW = imbinarize(imDataBG);
        figure;
        subplot(1,2,1); imshow(imDataBG_BW,[]);
        subplot(1,2,2); imshow(imDataBG,[]);
        
        %% Multiscale Closing Operation
        ws = round(5);
        se = strel('disk',ws);
        Io1 = imclose(imDataBG,se);
        
        ws = round(9); 
        se = strel('disk',ws);
        Io2 = imclose(imDataBG,se);
        
        Io = max(Io1, Io2);
        figure;
        subplot(1,3,1); imshow(Io1,[]);
        subplot(1,3,2); imshow(Io2,[]);
        subplot(1,3,3); imshow(Io, [])        
        
        
        %% Threshold Data - Version 2
        imDataBG_BW = imbinarize(Io);
        figure;
        subplot(1,3,1); imshow(imDataBG_BW,[]);
        subplot(1,3,2); imshow(imDataBG,[]);       
        
        %% 
        % Clear the Borders
        bw = imclearborder(imDataBG_BW);
        bw = imerode(bw,strel('disk',1));

        % Compute and 'Thin' Background Markers                
        distance_image = bwdist(bw);
        DL = watershed(distance_image, 4);
        basins = DL == 0;    
        
        gmag = imgradient(imDataBG);
        gmag2 = imimposemin(gmag, basins | imDataBG_BW);
        L = watershed(gmag2);
        L(~bw) = 0;
        
        rgb = label2rgb(L, 'jet');
        imshow(rgb); title('Watershed');
      
        
        %% Name and Save file
        filename = [subFolder(subfolderIdx).name '.pdf'];
        saveas(gcf,filename)
        
        % Retrieve Features and Normalize data according to pixelsize
        [feats,nFeats] = bwlabeln(bw);
        featsProp = regionprops(feats,bw,'ConvexArea','ConvexHull','Area','Centroid','MajorAxisLength','MinorAxisLength','Perimeter');
        
        %Pre-populate parameter
        lengths = [];
        area = [];
        convexarea = [];
        spreading = [];
        
        % Convert pixel to actual physical size
        for i = 1:1:length(featsProp)
            featsProp(i).MajorAxisLength = featsProp(i).MajorAxisLength*PhysicalSize;
            featsProp(i).Area = featsProp(i).Area*PhysicalSize^2;
            featsProp(i).ConvexArea = featsProp(i).ConvexArea*PhysicalSize^2;
            featsProp(i).Spreading = featsProp(i).ConvexArea/featsProp(i).Area;
            featsProp(i).Ratio = featsProp(i).Perimeter/featsProp(i).Area;
            spreading(i) = featsProp(i).Spreading;
            ratio(i) = featsProp(i).Ratio;
        end
        outputSpreading{folderIdx}{subfolderIdx-2} = spreading;
        meanSpreading{folderIdx}{subfolderIdx-2} = mean(outputSpreading{folderIdx}{subfolderIdx-2});
        medianSpreading{folderIdx}{subfolderIdx-2} = median(outputSpreading{folderIdx}{subfolderIdx-2});
        
        outputRatio{folderIdx}{subfolderIdx-2} = ratio;
        meanRatio{folderIdx}{subfolderIdx-2} = mean(outputRatio{folderIdx}{subfolderIdx-2});
        medianRatio{folderIdx}{subfolderIdx-2} = median(outputRatio{folderIdx}{subfolderIdx-2});
        numCells{folderIdx}{subfolderIdx-2} = length(featsProp);
        
        figure;
        hold on
        hist(spreading)
        filename = [subFolder(subfolderIdx).name '_spreading.pdf'];
        saveas(gcf,filename)
        close all
        
        save([subFolder(subfolderIdx).name 'stats.mat'],'featsProp');
        save([subFolder(subfolderIdx).name 'numCells.mat'],'numCells');
        
        clear featsProp;
    end
    numCell{folderIdx} = sum([numCells{folderIdx}{:}]);
end


%% Save all Output and their means
cd(savePath);
save('NumberCellsAll.mat','numCell');
save('Spreading.mat','outputSpreading');
save('meanSpreading.mat','meanSpreading');
save('medianSpreading.mat','medianSpreading');
save('Ratio.mat','outputRatio');
save('meanRatio.mat','meanRatio');
save('medianRatio.mat','medianRatio');
save('Conditions.mat','folderPath');

%% Save workspace
save([pathDate '.mat']);

disp('Output Stats Saved in');
cd(savePath)

%% Plot Data

figure;
hold on

y = num2cell(1:size(medianRatio,2));
x = cellfun(@(medianRatio, y) [y*ones(size(medianRatio,2),1) [medianRatio{:}]'], medianRatio, y, 'UniformOutput', false); 
X = vertcat(x{:});
boxplot(X(:,2), X(:,1))
title('Cell Spreading - Convex Area/Area');
hold off

savefig('Summary.fig')