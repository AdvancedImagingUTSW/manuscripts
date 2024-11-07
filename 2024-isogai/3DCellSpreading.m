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


%% Analysis of cell spreading in a 3D gel environment
% Prepare workspace by loading path and parameters
clc; clear;

%Name the folder to store final analysis
analysisFolder = 'stats summary';

%Define rawData path and list contents
dataDirectory = '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/raw/tisogai/fishSPIM/RPE/201120/';
ls(dataDirectory);

% Define the location for saving the analysis
savePath = '/project/bioinformatics/Danuser_lab/3Dmorphogenesis/analysis/tisogai/fishSPIM/RPE/201120_CK666/';
mkdirRobust(savePath)
cd(savePath);

%Make analysisFolder to organize data within savePath
mkdir(savePath,analysisFolder);
analysisSavePath = [savePath,analysisFolder];

%Determine directory contents of rawData
directoryContents = dir(dataDirectory);
for dirIdx = 3:1:length(directoryContents)
    folderPath{dirIdx-2} = directoryContents(dirIdx).name;
end

%Define pixel paramter and path based on imaging software
PhysicalSizeX = 0.406;
PhysicalSizeY = 0.406;  
PhysicalSizeZ = 1.060; 
pathDate = '';

%% Extract the Data
for folderIdx = 1:1:size(folderPath,2)
    disp(['Processing ' folderPath{folderIdx}]);
    folderPathLong = fullfile(dataDirectory,folderPath{folderIdx},pathDate);
    [numberCells] = numberOfStrings(folderPathLong,'Cell');
    
    for cellIdx = 1:1:numberCells
        disp(['Processing Cell Number ' num2str(cellIdx)]);
        [numberPositions] = numberOfStrings([folderPathLong filesep 'Cell' num2str(cellIdx)],'position');
        
        if(numberPositions>0)
            dataVector = 1:1:numberPositions;
        else
            dataVector = 1;
        end
        
        for positionIdx = dataVector
            filename = [folderPath{folderIdx} '_' num2str(cellIdx) '_' num2str(positionIdx) '.pdf'];
            isComplete = numberOfStrings(savePath,filename);
            if(isComplete==1)
                disp([filename ' complete']);
            elseif(isComplete==0)
                
                % Allocate Memory
                if(numberPositions>0)
                    disp(['Processing Position Number ' num2str(positionIdx)]);
                    imagePath = fullfile(folderPathLong,['Cell' num2str(cellIdx)],['position ' num2str(positionIdx)],['1_CH00_000000.tif']);
                else
                    imagePath = fullfile(folderPathLong,['Cell' num2str(cellIdx)],['1_CH00_000000.tif']);
                end
                imageInfo = imfinfo(imagePath);
                imData = zeros(imageInfo(1).Width,imageInfo(1).Height,length(imageInfo));
                
                %Load Data
                if(numberPositions>0)
                    imagePath = fullfile(folderPathLong,['Cell' num2str(cellIdx)],['position ' num2str(positionIdx)]);
                else
                    imagePath = fullfile(folderPathLong,['Cell' num2str(cellIdx)]);
                end
                 imData = load3DImage(imagePath,['1_CH00_000000.tif']);
                 imData = imresize3(imData,[size(imData,1)*PhysicalSizeX/PhysicalSizeZ size(imData,2)*PhysicalSizeX/PhysicalSizeZ round(size(imData,3))]);                %imData = load3DImage(imagePath,['1_CH00_' num2str(cellIdx) '.tif']);
                
                % Measure the Stats
                dataStat = [min(imData(:)) mean(imData(:)) std(imData(:)) max(imData(:))];
                
                % Threshold Data
                thresh = dataStat(2)+4*dataStat(3);
                binaryImage = imbinarize(imData,thresh);
                
                % Remove Small Connected Components
                % for non-spread cells, see a volume of 50,000, 44,000.  30k is
                % conservative.
                cellVolume = 400;
                binaryImage = bwareaopen(binaryImage,cellVolume);
                
                % Remove Border Components
                binaryImage = imclearborder(binaryImage);
                
                % Fill Holes
                labels = imfill(binaryImage,'holes');
                
                % Measure Region Properties
                stats = regionprops3(binaryImage,'BoundingBox','Centroid','ConvexVolume','PrincipalAxisLength','SurfaceArea','Volume');
                
                % Convert pixels into actual dimension in um
                for axis = 1:3
                    stats.PrincipalAxisLength(:,axis) = stats.PrincipalAxisLength(:,axis)*PhysicalSizeX;
                end
                
                % Save stats
                save([folderPath{folderIdx} '_' num2str(cellIdx) '_' num2str(positionIdx) '_stats.mat'],'stats');
                disp('Stats Saved');
                
                figure;
                hold on
                imshow(imcomplement(max(binaryImage,[],3)),[]);
                
                numberOfStats = size(stats);
                for statIdx = 1:1:numberOfStats(1);
                    regionLocation = stats.BoundingBox(statIdx,:);
                    lowerLeftCornerX = regionLocation(1);
                    lowerLeftCornerY = regionLocation(2);
                    rectangleWidth = regionLocation(4);
                    rectangleHeight = regionLocation(5);
                    rectangle('Position',[lowerLeftCornerX lowerLeftCornerY rectangleWidth rectangleHeight]);
                end
                filename = [folderPath{folderIdx} '_' num2str(cellIdx) '_' num2str(positionIdx) '.pdf'];
                imageName = [folderPath{folderIdx} '_' num2str(cellIdx) '_' num2str(positionIdx) '.tif'];
                stackWrite(binaryImage,imageName);
                saveas(gcf,filename)
                close all
            end
        end
    end
end
disp('Done');

%% Organize the Data
for folderIdx = 1:1:size(folderPath,2)
    finalData{folderIdx} = [];
    folderPathLong = fullfile(dataDirectory,folderPath{folderIdx},pathDate);
    [numberCells] = numberOfStrings(folderPathLong,'Cell');
    cd(savePath);
    for cellIdx = 1:1:numberCells
        [numberPositions] = numberOfStrings([folderPathLong filesep 'Cell' num2str(cellIdx)],'position');
        
        if(numberPositions>0)
            dataVector = 1:1:numberPositions;
        else
            dataVector = 1;
        end
        
        for positionIdx = dataVector
            disp(['Analysis Processing Folder # ' num2str(folderIdx) ' Cell # ' num2str(cellIdx) ' - Position # ' num2str(positionIdx)]);
            
            filename = [folderPath{folderIdx} '_' num2str(cellIdx) '_' num2str(positionIdx) '_stats.mat'];
            id = [folderPath{folderIdx}];
            load(filename);
            disp('Organizing stats......')
                statCell=cell(positionIdx,cellIdx);
                tmp=load([id '_' num2str(cellIdx) '_' num2str(positionIdx) '_stats.mat']);
                statCell{positionIdx}{cellIdx}=tmp.stats;
                statCell{positionIdx}{cellIdx}.elongation=(1-(statCell{positionIdx}{cellIdx}.PrincipalAxisLength(:,2)./statCell{positionIdx}{cellIdx}.PrincipalAxisLength(:,1)));
                data{positionIdx}{cellIdx} = [statCell{positionIdx}{cellIdx}.elongation];
                dataMean{positionIdx}{cellIdx} = mean(data{positionIdx}{cellIdx});
                dataMedian{positionIdx}{cellIdx} = median(data{positionIdx}{cellIdx});
                numCell{positionIdx}{cellIdx} = size(statCell{positionIdx}{cellIdx}.elongation,1);
        end
    end
    tempCell = [numCell{:,:}];
    allCell = sum(cell2mat(tempCell));
    cd(analysisSavePath);
    save([id 'data.mat'],'data');
    save([id 'dataMean.mat'],'dataMean');
    save([id 'dataMedian.mat'],'dataMedian');
    save([id 'numberCells.mat'],'allCell');
    display('Organized stats saved in')
    analysisSavePath
    clear data;
    clear dataMean;
    clear dataMedian;
    clear allCell;
    clear statCell;
end
save([pathDate '.mat']);

disp('Done');
