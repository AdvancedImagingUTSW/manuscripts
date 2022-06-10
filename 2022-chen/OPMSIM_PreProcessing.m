%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Pre-proceessing for OPM-SIM data
%           by Reto Fiolka, 09.09.2021
%
%           -Allow redo of cross-correlation matching
%           -Allow random starting timepoint, and estimate for that
%            timepoint correlation
%           -added automatic determination of z interpolation factor
%           -added axial shift correction via cross-correlation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clc
close all;
clear all;
cell=1;

folder='/archive/bioinformatics/Danuser_lab/Fiolka/Manuscripts/2021_OPMSIM/MatlabCode/ForNatureMethods/OPSIM_ExampleData/';

theta2=-60.35*pi/180;  % intial guess for relative rotation between dir 0 and dir 1
theta3=60.19*pi/180;
flipper=0;              % image flipping, needed for older data that had an addtional reflection in image path
zstep=0.32;
SM=round(zstep/0.16);


%% empirical parameters for image correlation
block=0;  %%block portions of image to be rotated
bsize=120;
bsizeX=120;

%% read in metadata
fileID = fopen([folder 'AcqInfo.txt'],'r');
Px=fgetl(fileID);
Py=fgetl(fileID);
Pz=fgetl(fileID);
fgetl(fileID);
fgetl(fileID);
fgetl(fileID);
timep=fgetl(fileID);
fclose(fileID);

xdim=str2num(Px(end-2:end));
ydim=str2num(Py(end-2:end));
if size(Pz,2)==14   %% for stacks of 100 planes or more
    nPlane=str2num(Pz(end-2:end));
else   % stack sizesfrom 10 to 99
    nPlane=str2num(Pz(end-1:end));
end
if size(timep,2)==15
    timepoints=num2str(timep(end-1:end));
else
    timepoints=num2str(timep(end));
end

timepoints=timepoints-1;



i=sqrt(-1);

X0=zeros(nPlane,ydim,xdim);
X1=zeros(nPlane,ydim,xdim);
X2=zeros(nPlane,ydim,xdim);

WF1=zeros(nPlane,512,xdim);
WF2=zeros(nPlane,512,xdim);

%%
count=1;
for n=0:1%timepoints
    
    
    if n<10
        name=['1_CH00_00000' num2str(n) '.tif'];
    else
        name=['1_CH00_0000' num2str(n) '.tif'];
    end
    
    if n>0
        clear X00 X01 X02 X10 X11 X12 X20 X21 X22
    end
    cor=0;
    
    %% read in data
    
    for k=1:nPlane
        if cor==0
            X00(k,:,:)=double(imread([folder name],k+0*nPlane-cor));
        else
            if k<cor+1
                X00(k,:,:)=double(imread([folder name],k+9*nPlane-cor));
            else
                X00(k,:,:)=double(imread([folder name],k+0*nPlane-cor));
            end
        end
        
        X01(k,:,:)=double(imread([folder name],k+nPlane+0*nPlane-cor));
        X02(k,:,:)=double(imread([folder name],k+2*nPlane+0*nPlane-cor));
    end
    
    for k=1:nPlane
        X10(k,:,:)=double(imread([folder name],k+3*nPlane-cor));
        X11(k,:,:)=double(imread([folder name],k+nPlane+3*nPlane-cor));
        X12(k,:,:)=double(imread([folder name],k+2*nPlane+3*nPlane-cor));
    end

    for k=1:nPlane
        X20(k,:,:)=double(imread([folder name],k+6*nPlane-cor));
        X21(k,:,:)=double(imread([folder name],k+nPlane+6*nPlane-cor)); 
        X22(k,:,:)=double(imread([folder name],k+2*nPlane+6*nPlane-cor));
    end
    
    
    %% deskew data
    
    X00=X00(:,5:end-5,5:end-5);
    X01=X01(:,5:end-5,5:end-5);
    X02=X02(:,5:end-5,5:end-5);
    X10=X10(:,5:end-5,5:end-5);
    X11=X11(:,5:end-5,5:end-5);
    X12=X12(:,5:end-5,5:end-5);
    X20=X20(:,5:end-5,5:end-5);
    X21=X21(:,5:end-5,5:end-5);
    X22=X22(:,5:end-5,5:end-5);
    
    if flipper
        X00=fliplr(X00);
        X01=fliplr(X01);
        X02=fliplr(X02);
        X10=fliplr(X10);
        X11=fliplr(X11);
        X12=fliplr(X12);
        X20=fliplr(X20);
        X21=fliplr(X21);
        X22=fliplr(X22);
    end
    
    X00=abs(permute(X00,[3,2,1]));
    X00=SupFun.shear3DinDim2(X00, 45, 1, zstep, 0.114, 0, 0);
    
    X01=abs(permute(X01,[3,2,1]));
    X01=SupFun.shear3DinDim2(X01, 45, 1, zstep, 0.114, 0, 0);
    
    X02=abs(permute(X02,[3,2,1]));
    X02=SupFun.shear3DinDim2(X02, 45, 1, zstep, 0.114, 0, 0);
    
    X10=abs(permute(X10,[3,2,1]));
    X10=SupFun.shear3DinDim2(X10, 45, 1, zstep, 0.114, 0, 0);
    
    X11=abs(permute(X11,[3,2,1]));
    X11=SupFun.shear3DinDim2(X11, 45, 1, zstep, 0.114, 0, 0);
    
    X12=abs(permute(X12,[3,2,1]));
    X12=SupFun.shear3DinDim2(X12, 45, 1, zstep, 0.114, 0, 0);
    
    X20=abs(permute(X20,[3,2,1]));
    X20=SupFun.shear3DinDim2(X20, 45, 1, zstep, 0.114, 0, 0);
    
    X21=abs(permute(X21,[3,2,1]));
    X21=SupFun.shear3DinDim2(X21, 45, 1, zstep, 0.114, 0, 0);
    
    X22=abs(permute(X22,[3,2,1]));
    X22=SupFun.shear3DinDim2(X22, 45, 1, zstep, 0.114, 0, 0);
    
    
    
    %% interpolation in z
    [xn,yn,zn]=size(X00);
    
    X00=imresize3(X00,[xn,yn,SM*zn],'linear');
    X01=imresize3(X01,[xn,yn,SM*zn],'linear');
    X02=imresize3(X02,[xn,yn,SM*zn],'linear');
    
    X10=imresize3(X10,[xn,yn,SM*zn],'linear');
    X11=imresize3(X11,[xn,yn,SM*zn],'linear');
    X12=imresize3(X12,[xn,yn,SM*zn],'linear');
    
    X20=imresize3(X20,[xn,yn,SM*zn],'linear');
    X21=imresize3(X21,[xn,yn,SM*zn],'linear');
    X22=imresize3(X22,[xn,yn,SM*zn],'linear');
    
    %% rotate data in coverslip frame
    theta=-45*pi/180;
    tform = affine3d([cos(theta) 0 sin(theta) 0; 0 1 0 0; -sin(theta) 0 cos(theta) 0; 0 0 0 1]);
    
    X00= imwarp(X00,tform,'linear');
    X01= imwarp(X01,tform,'linear');
    X02= imwarp(X02,tform,'linear');
    
    X10= imwarp(X10,tform,'linear');
    X11= imwarp(X11,tform,'linear');
    X12= imwarp(X12,tform,'linear');
    
    X20= imwarp(X20,tform,'linear');
    X21= imwarp(X21,tform,'linear');
    X22= imwarp(X22,tform,'linear');
    
    
    test=squeeze(max(X00,[],3));
    test=sum(test);
    ind=find(test>0);
    
    
    X00=X00(:,ind(1):ind(end),:);
    X01=X01(:,ind(1):ind(end),:);
    X02=X02(:,ind(1):ind(end),:);
    
    X10=X10(:,ind(1):ind(end),:);
    X11=X11(:,ind(1):ind(end),:);
    X12=X12(:,ind(1):ind(end),:);
    
    X20=X20(:,ind(1):ind(end),:);
    X21=X21(:,ind(1):ind(end),:);
    X22=X22(:,ind(1):ind(end),:);
    
    %%
    Y00=permute(X00,[1,3,2]);
    Y01=permute(X01,[1,3,2]);
    Y02=permute(X02,[1,3,2]);
    
    X10=permute(X10,[1,3,2]);
    X11=permute(X11,[1,3,2]);
    X12=permute(X12,[1,3,2]);
    
    Y20=permute(X20,[1,3,2]);
    Y21=permute(X21,[1,3,2]);
    Y22=permute(X22,[1,3,2]);
    
    %% crop center direction
    X10=X10(50:end-50,:,:);
    X11=X11(50:end-50,:,:);
    X12=X12(50:end-50,:,:);

    %% relative rotation of stacks to each other
    
    tform = affine3d([cos(theta3) -sin(theta3) 0 0; sin(theta3) cos(theta3) 0 0; 0 0 1 0; 0 0 0 1])
    X00= imwarp(Y00,tform,'linear');
    X01= imwarp(Y01,tform,'linear');
    X02= imwarp(Y02,tform,'linear');
    
    WF0=X00+X01+X02;
    WF1=X10+X11+X12;
    
    % cross correlation
    
    M1=squeeze(max(WF0,[],3));
    M2=squeeze(max(WF1,[],3));
    
    M1(1:bsize,:)=0;
    M1(end-bsize:end,:)=0;
    M1(:,1:(bsizeX-40))=0;
    M1(:,end-(bsizeX+60):end)=0;
    
    M2=M2(200:end-200,90:end-90); %Adjust carefully size of window; you can execute this section on its own again.
    M2=M2*100;
    M1=M1*100;
    if count==1
        
        c = normxcorr2(M2,M1);
        %figure; imshow(c,[])
        [ypeak, xpeak] = find(c==max(c(:)));
        
        % Compute translation from max location in correlation matrix
        yoffSet = ypeak-round(size(c,1)./2);
        xoffSet=xpeak-round(size(c,2)./2);
        max(c(:));
        
        
        x11=size(WF0,1);x12=size(WF0,2);
        x21=size(WF1,1);x22=size(WF1,2);
        
        off11=round((x11-x21)/2)+yoffSet;
        off21=round((x12-x22)/2)+xoffSet;
    end
    
    
    pad=400;
    X00=padarray(X00,[pad,pad]);
    X01=padarray(X01,[pad,pad]);
    X02=padarray(X02,[pad,pad]);
    WF0=X00+X01+X02;
    
    WF00=WF0(off11+pad:size(WF1,1)+off11-1+pad,off21+pad:size(WF1,2)+off21-1+pad,:);
    X00=X00(off11+pad:size(WF1,1)+off11-1+pad,off21+pad:size(WF1,2)+off21-1+pad,:);
    X01=X01(off11+pad:size(WF1,1)+off11-1+pad,off21+pad:size(WF1,2)+off21-1+pad,:);
    X02=X02(off11+pad:size(WF1,1)+off11-1+pad,off21+pad:size(WF1,2)+off21-1+pad,:);
    
    %figure;imshow(squeeze(max(WF00,[],3))-0.02*squeeze(max(WF1,[],3)),[])
    %%
    
    tform = affine3d([cos(theta2) -sin(theta2) 0 0; sin(theta2) cos(theta2) 0 0; 0 0 1 0; 0 0 0 1]);
    
    
    
    X20= imwarp(Y20,tform,'linear');
    X21= imwarp(Y21,tform,'linear');
    X22= imwarp(Y22,tform,'linear');
    WF2=X20+X21+X22;WF1=X10+X11+X12;
    
    if count==1
        M1=squeeze(max(WF2,[],3));
        M2=squeeze(max(WF1,[],3));
        M1(1:(bsize-50),:)=0;
        M1(end-bsize:end,:)=0;
        M1(:,1:bsizeX)=0;
        M1(:,end-(bsizeX+50):end)=0;
        
        
        M2=M2(150:end-250,60:end-60); % carefully adjust window; you can execute this section on its own.
        
        M2=M2*100;
        M1=M1*100;
        
        c = normxcorr2(M2,M1);
        %figure; imshow(c,[])
        [ypeak, xpeak] = find(c==max(c(:)));
        
        % Compute translation from max location in correlation matrix
        yoffSet = ypeak-round(size(c,1)./2)+50;
        xoffSet=xpeak-round(size(c,2)./2);
        max(c(:));
        
        
        x11=size(WF2,1);x12=size(WF2,2);
        x21=size(WF1,1);x22=size(WF1,2);
        
        off1=round((x11-x21)/2)+yoffSet;
        off2=round((x12-x22)/2)+xoffSet;
        
    end
    
    pad=400;
    X20=padarray(X20,[pad,pad]);
    X21=padarray(X21,[pad,pad]);
    X22=padarray(X22,[pad,pad]);
    WF2=X20+X21+X22;
    WF22=WF2(off1+pad:size(WF1,1)+off1-1+pad,off2+pad:size(WF1,2)+off2-1+pad,:);
    X20=X20(off1+pad:size(WF1,1)+off1-1+pad,off2+pad:size(WF1,2)+off2+pad-1,:);
    X21=X21(off1+pad:size(WF1,1)+off1-1+pad,off2+pad:size(WF1,2)+off2+pad-1,:);
    X22=X22(off1+pad:size(WF1,1)+off1-1+pad,off2+pad:size(WF1,2)+off2+pad-1,:);
    
    
    %% normalize stacks
    X00=X00/max(X00(:))*15000;
    X01=X01/max(X01(:))*15000;
    X02=X02/max(X02(:))*15000;
    
    X10=X10/max(X10(:))*15000;
    X11=X11/max(X11(:))*15000;
    X12=X12/max(X12(:))*15000;
    Y10=X10;
    Y11=X11;
    Y12=X12;
    
    X20=X20/max(X20(:))*15000;
    X21=X21/max(X21(:))*15000;
    X22=X22/max(X22(:))*15000;
    
    %%
    
    X10(:,:,1)=zeros(size(X10,1),size(X10,2));
    X10(:,:,2)=zeros(size(X10,1),size(X10,2));
    X10(:,:,3)=zeros(size(X10,1),size(X10,2));
    X10(:,:,4:end)=Y10(:,:,1:end-3);
    
    
    X11(:,:,1)=zeros(size(X10,1),size(X10,2));
    X11(:,:,2)=zeros(size(X10,1),size(X10,2));
    X11(:,:,3)=zeros(size(X10,1),size(X10,2));
    X11(:,:,4:end)=Y11(:,:,1:end-3);

    X12(:,:,1)=zeros(size(X10,1),size(X10,2));
    X12(:,:,2)=zeros(size(X10,1),size(X10,2));
    X12(:,:,3)=zeros(size(X10,1),size(X10,2));
    X12(:,:,4:end)=Y12(:,:,1:end-3);
    
    
    %%  axial correction
    
    r=squeeze(max(X01,[],2));
    g=squeeze(max(X11,[],2));
    rgbimg=cat(3,r/max(r(:)),g/max(g(:)),g/max(g(:)));
    figure;imshow(rgbimg,[]);
    
    c = normxcorr2(r,g);
    [ypeak, xpeak] = find(c==max(c(:)));
    yoffSet = ypeak-round(size(c,1)./2);
    xoffSet=xpeak-round(size(c,2)./2);
    X00=circshift(X00,xoffSet,3);
    X01=circshift(X01,xoffSet,3);
    X02=circshift(X02,xoffSet,3);
    
    %show axial mip of the two direction
    r=squeeze(max(X01,[],2));
    g=squeeze(max(X11,[],2));
    rgbimg=cat(3,r/max(r(:)),g/max(g(:)),g/max(g(:)));
    figure;imshow(rgbimg,[]);
    
    r=squeeze(max(X21,[],2));
    g=squeeze(max(X11,[],2));
    rgbimg=cat(3,r/max(r(:)),g/max(g(:)),g/max(g(:)));
    figure;imshow(rgbimg,[]);
    
    c = normxcorr2(r,g);
    [ypeak, xpeak] = find(c==max(c(:)));
    yoffSet = ypeak-round(size(c,1)./2);
    xoffSet=xpeak-round(size(c,2)./2);
    
    
    X20=circshift(X20,xoffSet,3);
    X21=circshift(X21,xoffSet,3);
    X22=circshift(X22,xoffSet,3);
    
    %show axial mip of the two direction
    r=squeeze(max(X21,[],2));
    g=squeeze(max(X11,[],2));
    rgbimg=cat(3,r/max(r(:)),g/max(g(:)),g/max(g(:)));
    figure;imshow(rgbimg,[]);
    
    %% cropping
    close all
    
    [xn,yn,zn]=size(X00);
    WF0=X00+X01+X02;
    WF1=X10+X11+X12;
    WF2=X20+X21+X22;
    
    if  count==1
        r=squeeze(max(WF0,[],3));
        g=squeeze(max(WF1,[],3));
        b=squeeze(max(WF2,[],3));
        rgbimg=cat(3,r/max(r(:)),g/max(g(:)),b/max(b(:)));
        figure('NumberTitle', 'off', 'Name', 'Please crop region of interest');
        [b rect]=imcrop(rgbimg);
    end
    
    
    for l=1:zn
        
        temp00=squeeze(X00(:,:,l));temp00=imcrop(temp00,rect);
        SIM00(:,:,l)=temp00;
        temp01=squeeze(X01(:,:,l));temp01=imcrop(temp01,rect);
        SIM01(:,:,l)=temp01;
        temp02=squeeze(X02(:,:,l));temp02=imcrop(temp02,rect);
        SIM02(:,:,l)=temp02;
        
        temp10=squeeze(X10(:,:,l));temp10=imcrop(temp10,rect);
        SIM10(:,:,l)=temp10;
        temp11=squeeze(X11(:,:,l));temp11=imcrop(temp11,rect);
        SIM11(:,:,l)=temp11;
        temp12=squeeze(X12(:,:,l));temp12=imcrop(temp12,rect);
        SIM12(:,:,l)=temp12;
        
        temp20=squeeze(X20(:,:,l));temp20=imcrop(temp20,rect);
        SIM20(:,:,l)=temp20;
        temp21=squeeze(X21(:,:,l));temp21=imcrop(temp21,rect);
        SIM21(:,:,l)=temp21;
        temp22=squeeze(X22(:,:,l));temp22=imcrop(temp22,rect);
        SIM22(:,:,l)=temp22;
        
    end
    
    %%
    if count==1
        WF0=X00+X01+X02;
        WF1=X10+X11+X12;
        WF2=X20+X21+X22;
        
        for l=1:zn
            
            W0(:,:,l)=imcrop(WF0(:,:,l),rect);
            W1(:,:,l)=imcrop(WF1(:,:,l),rect);
            W2(:,:,l)=imcrop(WF2(:,:,l),rect);
        end
        
        r=sum(W0,3);
        g=sum(W1,3);
        b=sum(W2,3);
        rgbimg=cat(3,r/max(r(:)),g/max(g(:)),b/max(b(:)));
        figure;imshow(rgbimg)
    end
    
    imgch0=SIM00+SIM01+SIM02;
    imgch1=SIM10+SIM11+SIM12;
    imgch2=SIM20+SIM21+SIM22;
    
    figure(1); imagesc(imgch0(:,:,round(size(imgch0,3)/2)));
    figure(2); imagesc(imgch1(:,:,round(size(imgch1,3)/2)));
    figure(3); imagesc(imgch2(:,:,round(size(imgch2,3)/2)));
    
    %%
    
    if count==1
        %% affine registration, multimodal correlation loss.
        downsample_factor = [16,8,4,2,1];
        mode = 'multimodal';
        iterations = [50,100,1000,1000,1000]; %[50,100,200,500,1000];
        type='affine';
        
        %% register 0 -> 1
        % step 1: affine
        [imgch0_reg,tform0] = SupFun.register3D_intensity_multiscale_affine(imgch1, imgch0, 1, 0, downsample_factor, mode, iterations, type);
        % step 2: non-rigid
        %[imgch0_reg,tform0_nonrigid] = register3D_demons(imgch1, imgch0_reg, [50, 100, 500], 1.);
        figure();
        imshowpair(squeeze(imgch1(:,:,round(size(imgch1,3)/2))), squeeze(imgch0_reg(:,:,round(size(imgch0_reg,3)/2))),'Scaling','joint');
        figure();
        imshowpair(squeeze(imgch1(:,round(size(imgch1,2)/2),:)), squeeze(imgch0_reg(:,round(size(imgch0_reg,2)/2),:)),'Scaling','joint');
        figure();
        imshowpair(squeeze(imgch1(round(size(imgch1,1)/2),:,:)), squeeze(imgch0_reg(round(size(imgch0_reg,1)/2),:,:)),'Scaling','joint');
        
        %% register 2->1
        [imgch2_reg,tform2] = SupFun.register3D_intensity_multiscale_affine(imgch1, imgch2, 1, 0, downsample_factor, mode, iterations, type);
        %[imgch2_reg,tform2_nonrigid] = register3D_demons(imgch1, imgch2_reg, [50, 100, 500], 1.);
        figure();
        imshowpair(squeeze(imgch1(:,:,round(size(imgch1,3)/2))), squeeze(imgch2_reg(:,:,round(size(imgch2_reg,3)/2))),'Scaling','joint');
        figure();
        imshowpair(squeeze(imgch1(:,round(size(imgch1,2)/2),:)), squeeze(imgch2_reg(:,round(size(imgch2_reg,2)/2),:)),'Scaling','joint');
        figure();
        imshowpair(squeeze(imgch1(round(size(imgch1,1)/2),:,:)), squeeze(imgch2_reg(round(size(imgch2_reg,1)/2),:,:)),'Scaling','joint');
        
        %%
        volCell={imgch0,imgch1,imgch2};
        
        
        %% Apply transform
        tformCell={affine3d(tform0'),affine3d(eye(4)),affine3d(tform2')};
        
    end
    volWF={imgch0,imgch1,imgch2};
    [stitchedVolume,volWF]=SupFun.stitchVolume(volWF,tformCell);
    
    volPH0={SIM00,SIM10,SIM20};
    [stitchedVolume,volPH0]=SupFun.stitchVolume(volPH0,tformCell);
    
    volPH1={SIM01,SIM11,SIM21};
    [stitchedVolume,volPH1]=SupFun.stitchVolume(volPH1,tformCell);
    
    volPH2={SIM02,SIM12,SIM22};
    [stitchedVolume,volPH2]=SupFun.stitchVolume(volPH2,tformCell);
    
    %% Saving volumes
    if count==1
        arrayfun(@(i) SupFun.stackWrite(uint16(volCell{i}),[folder 'reg_Widefield_' num2str(i-1) '-' num2str(n) '.tif']),1:numel(volCell));
    end
    arrayfun(@(i) SupFun.stackWrite(uint16(volPH0{i}),[folder 'reg_SIM' num2str(i-1) num2str(0) '-' num2str(n)  '.tif']),1:numel(volCell));
    arrayfun(@(i) SupFun.stackWrite(uint16(volPH1{i}),[folder 'reg_SIM' num2str(i-1) num2str(1) '-' num2str(n) '.tif']),1:numel(volCell));
    arrayfun(@(i) SupFun.stackWrite(uint16(volPH2{i}),[folder 'reg_SIM' num2str(i-1) num2str(2) '-' num2str(n) '.tif']),1:numel(volCell));
    disp('done saving data')
    count=count+1;
    
end
