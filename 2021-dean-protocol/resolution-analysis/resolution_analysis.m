%% Resolution Analysis
plot_figure = 1;
export_image = 1;
export_CSV = 1;
export_line_profiles = 1;

data_path = '/archive/MIL/dean/20210623_code_optimization_for_git';
image_name = '1_CH00_000000.tif';

imPath = fullfile(data_path,image_name);
lateral_pixel_size = 115;
axial_pixel_size = 115;
image_threshold = 350/2^16;

cd(data_path);
export_path = fullfile(data_path,'ROIs');
if(~isdir(export_path))
    mkdir(export_path);
end

%% 
% Get the information regarding the images we are going to work with.
image_info = imfinfo(imPath);

% Load the Data
number_Z = size(image_info,1);
for z_Idx = 1:1:number_Z
    im_data(:,:,z_Idx) = imread(imPath,z_Idx);
end
disp('Data Loaded');

% Threshold the Data.  Here we do something stupid and just do an
% arbitrary hard threshold.  These values are defined for each
% dataset above.
im_BW = imbinarize(im_data, image_threshold);

% Get the information from the regions that were above the
% threshold.
stats = regionprops(im_BW,'centroid');
disp('Particles Identified')

%% Go through each bead and measure its FWHM
% Initialize the finalResults structure.
number_beads = size(stats,1);
final_results = struct('fwhmX',[],'fwhmY',[],'fwhmZ',[],'xR2',[],'yR2',[],'zR2',[]);


% Iterate through each bead that was identified.
for bead_idx = 1:1:number_beads
    
    % Display Progress
    disp([num2str(bead_idx) ' out of ' num2str(number_beads) ' evaluated']);
    
    % Identify bead centroid around round to nearest integer value.
    % The roiRadius specifies how big of an ROI we are going to
    % take.
    roi_radius = 7;
    start_idx = [];
    end_idx = [];
    start_idx(1) = stats(bead_idx).Centroid(2)-roi_radius;
    end_idx(1) = stats(bead_idx).Centroid(2)+roi_radius;
    
    start_idx(2) = stats(bead_idx).Centroid(1)-roi_radius;
    end_idx(2) = stats(bead_idx).Centroid(1)+roi_radius;
    
    start_idx(3) = stats(bead_idx).Centroid(3)-roi_radius;
    end_idx(3) = stats(bead_idx).Centroid(3)+roi_radius;
    
    start_idx = round(start_idx);
    end_idx = round(end_idx);
    
    % Make sure particle is not touching the edge.
    % If negative number, too small.
    % If the quotient is greater than 1, too big.
    if(all(start_idx > 0) && all(end_idx./size(im_BW) < 1))
        
        % Create local ROI that encompasses the bead.
        im_ROI = im_data(start_idx(1):end_idx(1),start_idx(2):end_idx(2),start_idx(3):end_idx(3));
        
        % Fine the brighest voxel, and then take line profiles
        % through it from each dimension.
        [~, maxIdx]=max(im_ROI(:));
        [ypeak,xpeak,zpeak]= ind2sub(size(im_ROI),maxIdx);
        lineX = double(squeeze(im_ROI(ypeak,:,zpeak)));
        lineY = double(squeeze(im_ROI(:,xpeak,zpeak)));
        lineZ = double(squeeze(im_ROI(ypeak,xpeak,:)));
        
        % Background subtract and normalize.
        lineZ=lineZ-min(lineZ); lineZ=lineZ/max(lineZ);
        lineX=lineX-min(lineX); lineX=lineX/max(lineX);
        lineY=lineY-min(lineY); lineY=lineY/max(lineY);
        
        % Specify Fit Options: f(x) =  a1*exp(-((x-b1)/c1)^2)
        fo = fitoptions('gauss1');
        fo.Lower=[0, 0.25*size(lineX,2), 0];
        fo.Upper=[1,0.75*size(lineX,2), 20];
        
        % Fit the Data
        x = 1:size(lineX,2);  x = x';
        [x_results, xgof] = fit(x,lineX','gauss1',fo);
        y = 1:size(lineY,1);  y = y';
        [y_results, ygof] = fit(y,lineY,'gauss1',fo);
        z = 1:size(lineZ,1);  z = z';
        [z_results, zgof] = fit(z,lineZ,'gauss1',fo);
        
        % Convert to FWHM
        scale_factor = 2*sqrt(2*log(2));
        fwhm_X = x_results.c1*lateral_pixel_size*scale_factor/sqrt(2);
        fwhm_Y = y_results.c1*lateral_pixel_size*scale_factor/sqrt(2);
        fwhm_Z = z_results.c1*axial_pixel_size*scale_factor/sqrt(2);
        
        % Place the FWHM results into a structure
        final_results.fwhmX(bead_idx)= fwhm_X;
        final_results.fwhmY(bead_idx)= fwhm_Y;
        final_results.fwhmZ(bead_idx)= fwhm_Z;
        
        % Place the R^2  results into a structure
        final_results.xR2(bead_idx)=xgof.rsquare;
        final_results.yR2(bead_idx)=ygof.rsquare;
        final_results.zR2(bead_idx)=zgof.rsquare;
        
        % Place the Location results into a structure
        final_results.x(bead_idx)=stats(bead_idx).Centroid(2);
        final_results.y(bead_idx)=stats(bead_idx).Centroid(1);
        final_results.z(bead_idx)=stats(bead_idx).Centroid(3);
        
        % If you want to export each bead as a sub-volume.  Useful for
        % loading each bead into ImageJ or an equivalent program.
        if(export_image == 1)
            quality_of_fit = [final_results.xR2(bead_idx) final_results.yR2(bead_idx) final_results.zR2(bead_idx)];
            
            % Only export beads with a good quality of fit.
            if(all(quality_of_fit > 0.95))
                export_path = fullfile(data_path,'ROIs');
                mkdir(export_path);
                filename = [num2str(fwhm_X) '_' num2str(fwhm_Y) '_' num2str(fwhm_Z) '.tif'];
                write3Dtiff(im_ROI, fullfile(export_path,filename));
            end
        end
        
        % Plot the Results.  Useful for visual inspection on a PDF.
        if(plot_figure == 1)
            quality_of_fit = [final_results.xR2(bead_idx) final_results.yR2(bead_idx) final_results.zR2(bead_idx)];
            if(all(quality_of_fit > 0.95))
                filename = [num2str(fwhm_X) '_' num2str(fwhm_Y) '_' num2str(fwhm_Z) '.pdf'];
                
                h = figure;
                subplot(3,1,1);
                plot(x_results,x,lineX);
                title('X-Axis');
                hLeg = legend('example');
                set(hLeg,'visible','off');
                
                subplot(3,1,2);
                plot(y_results,y,lineY);
                title('Y-Axis');
                hLeg = legend('example');
                set(hLeg,'visible','off');
                
                subplot(3,1,3);
                plot(z_results,z,lineZ);
                title('Z-Axis');
                hLeg = legend('example');
                set(hLeg,'visible','off');
                
                saveas(gcf, fullfile(export_path,filename));
                close all
            end
        end
        
        % Export a CSV file of the data.  Useful for double-checking curve
        % fitting.
        if(export_CSV == 1)
            quality_of_fit = [final_results.xR2(bead_idx) final_results.yR2(bead_idx) final_results.zR2(bead_idx)];
            if(all(quality_of_fit > 0.95))
                % X-Axis
                output = [0:1:size(lineX,2)-1];
                output(2,:) = lineX;
                output = output';
                filename = [num2str(fwhm_X) '_' num2str(fwhm_Y) '_' num2str(fwhm_Z) '_X.txt'];
                writematrix(output,fullfile(export_path,filename));
                clear output
                
                % Y-Axis
                output = [0:1:size(lineY,1)-1];
                output(2,:) = lineY;
                output = output';
                filename = [num2str(fwhm_X) '_' num2str(fwhm_Y) '_' num2str(fwhm_Z) '_Y.txt'];
                writematrix(output,fullfile(export_path,filename));
                clear output
                
                % Z-Axis
                output = [0:1:size(lineZ,1)-1];
                output(2,:) = lineZ;
                output = output';
                filename = [num2str(fwhm_X) '_' num2str(fwhm_Y) '_' num2str(fwhm_Z) '_Z.txt'];
                writematrix(output,fullfile(export_path,filename));
                clear output
            end
        end
        
        % Export a 1D tiff file.  Again, useful for evaluating curve fit.
        if(export_line_profiles == 1)
            quality_of_fit = [final_results.xR2(bead_idx) final_results.yR2(bead_idx) final_results.zR2(bead_idx)];
            if(all(quality_of_fit > 0.95))
                filename = [num2str(fwhm_X) '_' num2str(fwhm_Y) '_' num2str(fwhm_Z) '_Z.tiff'];
                imwrite(lineZ, fullfile(export_path,filename));
                
                filename = [num2str(fwhm_X) '_' num2str(fwhm_Y) '_' num2str(fwhm_Z) '_X.tiff'];
                imwrite(lineX, fullfile(export_path,filename));
                
                filename = [num2str(fwhm_X) '_' num2str(fwhm_Y) '_' num2str(fwhm_Z) '_Y.tiff'];
                imwrite(lineY, fullfile(export_path,filename));
            end
        end
        % Clear variables to make sure that they aren't carried through the
        % for loop.
        clear x_results xgof y_results ygof z_results zgof imROI lineX lineY lineZ fwhmX fwhmY fwhmZ
    end
end

