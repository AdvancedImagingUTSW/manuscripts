%% small program to simulate projections of a rotated and sheared image.

close all

theta = -45;   %rotation and shearing angle

Im=double(imread('cameraman.tif'));

%rotate the image
Im2=imrotate(Im,theta);

[x,y]=size(Im2);

figure; imshow(Im2,[])

%project the rotated image
line=sum(Im2,1);


%compute matrix for affine transformation
tform = affine2d([ ...
    1 0 0;...
    tand(theta) 1 0; ...
    0 0 1]);

%shear the image
Im3=imwarp(Im,tform);

[x2,y2]=size(Im3);

figure;imshow(Im3,[])

%project the sheared image
line2=sum(Im3,1);
XI=[1:length(line2)];
XQ=linspace(1,length(line2),y);
%stretch the projected image to same size as line profile from rotated and
%projected image
line3 = interp1(XI,line2,XQ); 

%normalize line profile intensities
line=line/max(line);line3=line3/max(line3);
figure;plot(line3);hold on; plot(line,'r-');legend('shear project', 'rotate project')


