%% Code to get the curved projection image. 
% The image after de-scanning won't move, so we can control the pixel
% numbers and shear speed. But there can be a base speed according to the
% scanning range.
clc; 
clear;
close all;

% read calculated curve waveform
imagePath = 'F:\Bingying\Curvedprops_mode\embyro1';
wavefrontName='scanning_curve_0901-1.mat';
infoName = 'info_projection.txt';
filepath=fullfile(imagePath,wavefrontName);
load(filepath)

delay = 0; % percentage of the period. Add delay to scanning and shearing. When it's negative, add it to camera trigger.
delay2 = 0;

time_interval = 9.5; %unit s. Add time interval for time lapses. Acquisition time included. 
timepoint= 200; %if timepoint= 0, continuous mode

shear_factor= 1.41;  %(1+tand(45)*tand(viewing_angle)); 1, descrewed view; 2, around top-down view;
shear_pixel2voltage= 7.9813e-4;% voltage per pixel
xgalvo_scan2voltage= 1/100; %voltage per 100um


direction= 1; % /__/ left to right is 1. for opmv2, from negtive value to positive value
plot_ouput= 1; % 1, plot ouput control signal 

Vrs= 20; % Rolling shutter

recommand_Ypixelnumber = length(scanning_curve);

disp(['------------------------------------------------------------------------']);
disp(['Please set the pixel number in Y as ' num2str(recommand_Ypixelnumber) ', and Y0 as ' num2str((2048-recommand_Ypixelnumber)/2) '.']);

% Setting limitation for HCImageLive

Vn= 1514; %recommand_Ypixelnumber; % number_vertical_pixels
disp(['------------------------------------------------------------------------']);
disp(['When Y is ' num2str(Vn) ', Y0 should be ' num2str((2048-Vn)/2) '.']);

Hmin = 9.74436E-6; % minimum value, 9.74436E-6,unit seconds
Emin = (Vrs+Vn+10)*Hmin*1000; % minimum exposure time, unit ms
x_galvo_duty_cycle = 0.95;
FRmax = 1/Emin*1000*x_galvo_duty_cycle; % fatest frame rateHz

disp(['Maximum frame rate is ' num2str(FRmax) ' Hz.']);
disp(['Please set the frame rate less than that.']);

% Set line interval HCImageLive

imaging_frequency = 10; % must be smaller than FRmax
H = 1/imaging_frequency/(Vrs+Vn+10)*x_galvo_duty_cycle;
exposure= Vrs*H;
disp(['------------------------------------------------------------------------']);
disp(['Set the Line Interval as ' num2str(H*10^6) ' us']);
disp(['Set the exposure as ' num2str(exposure*10^3) ' ms']); % might not be necessary

sampling_rate = 100000; % Hz. limit to 5000 for USB-DAQ. 1 MHz for NI PCIe-6738
%% Initialize the Acquisition

% Configure the Camera

% Standard Scan: Camera Link - Orca Flash 4.0 v3
% External trigger mode (Edge trigger / Level trigger)
assert( Vn <= 2048); % number_vertical_pixels
imaging_duration= timepoint/imaging_frequency;
disp(['------------------------------------------------------------------------']);
disp(['The Imaging Duration is ' num2str(imaging_duration) ' s']);
disp(['The Rolling Shutter width is ' num2str(Vrs) ' pixels']);
disp(['The specified Imaging Frequency is ' num2str(imaging_frequency) ' Hz']);
disp(['A total of ' num2str(timepoint) ' Images will be Captured']);

% Load DAQ
daq_list = daqlist("ni");

%   "PXI6259"    "National Instruments(TM) PXI-6259"    "PXI-6259"    [1×1 daq.ni.PXIModule] Instruments(TM) PXI-6733"    "PXI-6733"    [1×1 daq.ni.PXIModule]
%   16-Bit, 8-Channel, 1 MS/s PXI Analog Output Module
%   [x_galvo' shear_galvo' laserShutter' edge_trigger'];
%   [ao0      ao1          Port0/Line1   Port0/Line0]; 
    daq_object = daq("ni");
    addoutput(daq_object,"Dev5","ao0","Voltage");
    %ch.TerminalConfig = 'On Demand';
    addoutput(daq_object,"Dev5","ao1","Voltage");
    %ch.TerminalConfig = 'On Demand';
    addoutput(daq_object,"Dev5",'Port0/Line1','Digital');
    %ch.TerminalConfig = 'On Demand';
    daq_object.Rate = sampling_rate;
    addoutput(daq_object,"Dev5",'Port0/Line0','Digital');
    
    
%% Generate Waveforms
% Time
t = 0:(1/sampling_rate):1/imaging_frequency;
camera_delay= ceil(10*H/(1/sampling_rate)); % delay and jitter 10H

%Add delay
delay = round(delay*length(t)/100);
delay2 = round(delay2*length(t)/100);

% Camera Trigger
ttl_amplitude = 1; %trigger voltage

edge_trigger = zeros(size(t));
edge_trigger(1:20) = ttl_amplitude;


%x_galvo, resample x-galvo curve and add the flyback, 
l_=1:length(scanning_curve);
l=linspace(1,length(scanning_curve),ceil(length(t)*x_galvo_duty_cycle));
XX=interp1(l_,scanning_curve,l,'spline');  % use whatever method suits you the best

flyback= linspace(XX(end),XX(1),length(t)-length(XX));

x_galvo= [XX,flyback];
x_galvo= x_galvo*xgalvo_scan2voltage;


% shear_galvo, resample to get shear_curve_shift same data size
shear_galvo_base =0.5*(sawtooth(2*pi*imaging_frequency*t,x_galvo_duty_cycle)); 
shear_galvo_base =shear_galvo_base*(Vn+Vrs);

YY=interp1(l_,shearing_curve,l,'spline');  % use whatever method suits you the best
shear_flyback= linspace(YY(end),YY(1),length(t)-length(YY));

shear_curve_shift= [YY, shear_flyback];

shear_galvo= (shear_galvo_base+shear_curve_shift)*shear_pixel2voltage;

shear_galvo_base= shear_galvo_base*shear_pixel2voltage;
shear_curve_shift=shear_curve_shift*shear_pixel2voltage;

%laser shutter

laserShutter=zeros(size(x_galvo));
shutter_end= round(length(laserShutter)*x_galvo_duty_cycle);
laserShutter(camera_delay:shutter_end+delay)=1;


if delay>0
    shear_galvo= [repmat(shear_galvo(1), 1, delay),shear_galvo];
    shear_galvo= shear_galvo(1:length(t));
    
    shear_galvo_base= [repmat(shear_galvo_base(1), 1, delay),shear_galvo_base];
    shear_galvo_base= shear_galvo_base(1:length(t));

    x_galvo= [repmat(x_galvo(1), 1, delay),x_galvo];
    x_galvo= x_galvo(1:length(t));
  
    if delay2>0
        shear_galvo= [repmat(shear_galvo(1), 1, delay2),shear_galvo];
        shear_galvo= shear_galvo(1:length(t));
    else
        if delay2<0
            x_galvo= [repmat(x_galvo(1), 1, -delay2),x_galvo];
            x_galvo= x_galvo(1:length(t));
        end
    end

        
else
    if delay<0     
        edge_trigger= [repmat(edge_trigger(end), 1, -delay),edge_trigger];
        edge_trigger= edge_trigger(1:length(t));
    end
end
  

%% for timepoint >1
if timepoint >=1
    if time_interval ==0
        x_galvo = repmat(x_galvo,[1,timepoint]);
        shear_galvo = repmat(shear_galvo,[1,timepoint]);
        edge_trigger=repmat(edge_trigger,[1,timepoint]);
        laserShutter=zeros(size(x_galvo));
        laserShutter(camera_delay:end-camera_delay)=1;
        
        output_matrix = [x_galvo' shear_galvo' laserShutter' edge_trigger'];
        write(daq_object,output_matrix);
        
    else
        for i= 1:timepoint
            tic
            output_matrix = [x_galvo' shear_galvo' laserShutter' edge_trigger'];
            write(daq_object,output_matrix);  
            pause(time_interval)
            toc
        end
    end
end
    

%% plot the output

if plot_ouput
    plot(x_galvo);
    hold on;
    plot(shear_galvo)
    plot(edge_trigger)
    plot(laserShutter)
    plot(shear_galvo_base)
    hold off;
end

leng=length(x_galvo);


%% for continuous

if timepoint ==0
    hWaitbar = waitbar(0, 'Iteration 1', 'Name', 'Solving problem','CreateCancelBtn','delete(gcbf)');
    shear_galvo_base=shear_galvo_base(1:leng);
    for i=1:300
        % Some long taking computation
        output_matrix = [x_galvo' shear_galvo' laserShutter' edge_trigger'];
        write(daq_object,output_matrix);
        % Check
        drawnow;

        if ~ishandle(hWaitbar)
            % Stop the if cancel button was pressed
            disp('Stopped by user');
            break;
        else
            % Update the wait bar
            waitbar(i/5,hWaitbar, ['Iteration ' num2str(i)]);
        end
    end
end
  
write(daq_object,[0 0 0 0]);


%% save imaging info

if timepoint ~=0
    filepath= fullfile(imagePath,infoName);
    str1= ['x Galvo scan range: ' num2str(scanrange(1)) ' um to ' num2str(scanrange(2)) ' um.'];
    str2= [str1 newline 'Rolling shutter width is ' num2str(Vrs) ' pixels.'];
    str3= [str2 newline 'Shear factor is ' num2str(shear_factor) '.'];

    writematrix(str2, filepath);
end
