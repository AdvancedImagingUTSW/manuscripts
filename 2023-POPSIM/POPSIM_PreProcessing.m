 clc; clear all; close all

        
        
        
folder='/archive/bioinformatics/Danuser_lab/Fiolka/Manuscripts/POPSIM/CodeToShare/ExampleData/'

%% cross-correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time=[1];                % timepoint to be processed
cro=0;                    % allow manual cropping with crop=1
rector=[64 125 652 578];  % corner points of fixed cropping area
affine=1;                 % use cross correlation for coarse estimation with affine=1; affine=0 uses hardcoded values

% intial guess for relative rotation angles
theta2=-59;  
theta3=60.5;


    
tic
timepoint= time;
  for n=time(1):time(end)
     
      if n<10  
      name=['1_CH00_00000' num2str(n) '.tif'];
      else 
      name=['1_CH00_0000' num2str(n) '.tif']
      end
          
  for k=1:9
        
temp00(k,:,:)=(imread([folder name],k));

  end
[x1,x2]=size(squeeze(temp00(1,:,:)));

 for k=1:9
 tempA=zeros(x1,x2);
 temp=squeeze(temp00(k,:,:));

temp=imresize(temp,[floor(x1./sqrt(2)),x2]);
[x11,x22]=size(temp);

diff=floor((x1-x11)/2);

tempA(diff+1:diff+x11,:)=temp;
% 
 FinalImage2(k,:,:)=(tempA);
% 
 end;

X00=squeeze(FinalImage2(1,:,:));
X01=squeeze(FinalImage2(2,:,:));
X02=squeeze(FinalImage2(3,:,:));

X10=squeeze(FinalImage2(4,:,:));
X11=squeeze(FinalImage2(5,:,:));
X12=squeeze(FinalImage2(6,:,:));

X20=squeeze(FinalImage2(7,:,:));
X21=squeeze(FinalImage2(8,:,:));
X22=squeeze(FinalImage2(9,:,:));




bl=size(X00);
blub=sort(bl,'ascend');
sdim=blub(1);

X00=X00(1:sdim,1:sdim);
X01=X01(1:sdim,1:sdim);
X02=X02(1:sdim,1:sdim);

X10=X10(1:sdim,1:sdim);
X11=X11(1:sdim,1:sdim);
X12=X12(1:sdim,1:sdim);

X20=X20(1:sdim,1:sdim);
X21=X21(1:sdim,1:sdim);
X22=X22(1:sdim,1:sdim);


WF0= X00+X01+X02;
WF1= X10+X11+X12;
WF2= X20+X21+X22;

if n == time(1)
     
    if affine==1
    fun=@(x)SupFun.AffImFcn(x,WF0,WF1);
    a0=[1,1,theta2];
    x = fminsearch(fun,a0)
       if theta2-2.5<x(3) & x(3)<theta2+2.5
           ftheta2 = x(3)
       else
           ftheta2 = theta2
       end
   
    
    
    fun=@(x)SupFun.AffImFcn(x,WF2,WF1);
    a0=[1,1,theta3];
    x = fminsearch(fun,a0)
       if theta3-2.5<x(3) & x(3)<theta3+2.5
           ftheta3 = x(3)
       else
           ftheta3 = theta3
       end
    else
    
        ftheta2 = theta2
        ftheta3 = theta3
    end
end

X00= imrotate(X00,theta2,'bicubic','crop');
X01= imrotate(X01,theta2,'bicubic','crop');
X02= imrotate(X02,theta2,'bicubic','crop');

X20= imrotate(X20,theta3,'bicubic','crop');
X21= imrotate(X21,theta3,'bicubic','crop');
X22= imrotate(X22,theta3,'bicubic','crop');

WF0= X00+X01+X02;
WF1= X10+X11+X12;
WF2= X20+X21+X22;

%% cross correlation
bsize=100;bsizeX=100
M1=WF0;
M2=WF1;
 M1(1:bsize,:)=0;
 M1(end-bsize:end,:)=0;
 M1(:,1:(bsizeX-40))=0;
 M1(:,end-(bsizeX+60):end)=0;

M2=M2(200:end-200,90:end-90); %Adjust carefully size of window; you can execute this section on its own again.
M2=M2*100;
M1=M1*100;


    c = normxcorr2(M2,M1);
    %figure; imshow(c,[])
max(c(:))
    
        [ypeak, xpeak] = find(c==max(c(:)));
    % Compute translation from max location in correlation matrix
    yoffSet = ypeak-round(size(c,1)./2)
    xoffSet=xpeak-round(size(c,2)./2)
    max(c(:))
    
    
x11=size(WF0,1);x12=size(WF0,2);
x21=size(WF1,1);x22=size(WF1,2);

off11=round((x11-x21)/2)+yoffSet
off21=round((x12-x22)/2)+xoffSet
if affine==0
 off11= -20; %smaller value: blue down
 off21=10;  %small0r value: blue to the right
end


WF00=circshift(WF0,[-off11,-off21]);
X00=circshift(X00,[-off11,-off21]);
X01=circshift(X01,[-off11,-off21]);
X02=circshift(X02,[-off11,-off21]);

    g=zeros(size(WF1));
    r=squeeze(max(WF1,[],3));
    b=squeeze(max(WF00,[],3));
    rgbimg=cat(3,r/max(r(:)),g/max(g(:)),b/max(b(:)));




WF2=X20+X21+X22;
WF1=X10+X11+X12;


M1=squeeze(max(WF2,[],3));
M2=squeeze(max(WF1,[],3));
M1(1:(bsize-50),:)=0;
M1(end-bsize:end,:)=0;
M1(:,1:bsizeX)=0;
M1(:,end-(bsizeX+50):end)=0;


M2=M2(50:end-50,50:end-50); % carefully adjust window; you can execute this section on its own.

M2=M2*100;
M1=M1*100;
% M2(1:100,:)=0;
% M2(end-250:end,:)=0;

    c = normxcorr2(M2,M1);
    %max(c(:))
    %figure; imshow(c,[])
        [ypeak, xpeak] = find(c==max(c(:)));
    % Compute translation from max location in correlation matrix
    yoffSet = ypeak-round(size(c,1)./2);
    xoffSet=xpeak-round(size(c,2)./2)
  
    
    
x11=size(WF2,1);x12=size(WF2,2);
x21=size(WF1,1);x22=size(WF1,2);

off1=round((x11-x21)/2)+yoffSet
off2=round((x12-x22)/2)+xoffSet
if affine==0
 off1=-20;%smaller value: blue down
 off2=-15;%smaller value: blue to the right
end

WF22=circshift(WF2,[-off1,-off2]);
X20=circshift(X20,[-off1,-off2]);
X21=circshift(X21,[-off1,-off2]);
X22=circshift(X22,[-off1,-off2]);
    r=squeeze(max(WF1,[],3));
    g=zeros(size(WF1));
    b=squeeze(max(WF22,[],3));
    rgbimg=cat(3,r/max(r(:)),g/max(g(:)),b/max(b(:)));
 
 
    r=squeeze(max(WF00,[],3));
    g=squeeze(max(WF1,[],3));
    b=squeeze(max(WF22,[],3));
    rgbimg=cat(3,r/max(r(:)),g/max(g(:)),b/max(b(:)));
    

     
  toc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 
%% readfiles

    tic

SIM00=uint16(X00);
SIM01=uint16(X01);
SIM02=uint16(X02);

SIM10=uint16(X10);
SIM11=uint16(X11);
SIM12=uint16(X12);

SIM20=uint16(X20);
SIM21=uint16(X21);
SIM22=uint16(X22);

imgch0=uint16(SIM00+SIM01+SIM02);
imgch1=uint16(SIM10+SIM11+SIM12);
imgch2=uint16(SIM20+SIM21+SIM22);



%% cropping 
if n==timepoint(1)
    r=imgch0;
    g=imgch1;
    b=imgch2;
    rgbimg=cat(3,imadjust(r),imadjust(g),imadjust(b));
    figure(4); 
    if cro
       
    [b rect]=imcrop(rgbimg);
    rect
    else
        rect=rector;
    end
end
    SIM00= imcrop(SIM00,rect);
    SIM01= imcrop(SIM01,rect);
    SIM02= imcrop(SIM02,rect);
    
    SIM10= imcrop(SIM10,rect);
    SIM11= imcrop(SIM11,rect);
    SIM12= imcrop(SIM12,rect);
    
    SIM20= imcrop(SIM20,rect);
    SIM21= imcrop(SIM21,rect);
    SIM22= imcrop(SIM22,rect);
    
    imgch0=SIM00+SIM01+SIM02;
    imgch1=SIM10+SIM11+SIM12;
    imgch2=SIM20+SIM21+SIM22;  
    imgchref=imgch0+imgch1+imgch2;
    
% figure(1); imagesc(max(imgch0,[],3))
% figure(2); imagesc(max(imgch1,[],3))
% figure(3); imagesc(max(imgch2,[],3))
% close(figure(4))
if n==timepoint(1)
%% affine registration, multimodal correlation loss. 
downsample_factor = [16,8,4,2,1]; % translation - [4,2,1] -> 10000 each, affine - [1] -> 10000
mode = 'multimodal'; 
iterations = [50,100,1000,1000,1000]; % originally [50,100,200,500,1000]
type='affine'; 

%% register 0 -> 1
% step 1: affine 
[imgch0_reg,tform0] = SupFun.register2D_intensity_multiscale_affine(imgch1, imgch0, 1, 0, downsample_factor, mode, iterations, type);
% step 2: non-rigid
%[imgch0_reg,tform0_nonrigid] = register3D_demons(imgch1, imgch0_reg, [50, 100, 500], 1.);
% figure(1)
% imshowpair(imgch1, imgch0_reg)

%% register 2 -> 1
% step 1: affine 
[imgch2_reg,tform2] = SupFun.register2D_intensity_multiscale_affine(imgch1, imgch2, 1, 0, downsample_factor, mode, iterations, type);
% step 2: non-rigid
%[imgch0_reg,tform0_nonrigid] = register3D_demons(imgch1, imgch0_reg, [50, 100, 500], 1.);
%  figure(2)
% imshowpair(imgch1, imgch2_reg)

    r= imgch1;
    g= imgch0_reg;
    b= imgch2_reg;
    rgbimg=cat(3,imadjust(r),imadjust(g),imadjust(b));
%     figure(3);imagesc(rgbimg)
 
    Rfixed = imref2d(size(imgch1));   
    

end
%% apply transformation

%affine3d(tform2')
r_SIM00 = imwarp(SIM00,affine2d(tform0'),'OutputView',Rfixed);
r_SIM01 = imwarp(SIM01,affine2d(tform0'),'OutputView',Rfixed);
r_SIM02 = imwarp(SIM02,affine2d(tform0'),'OutputView',Rfixed);

r_SIM20 = imwarp(SIM20,affine2d(tform2'),'OutputView',Rfixed);
r_SIM21 = imwarp(SIM21,affine2d(tform2'),'OutputView',Rfixed);
r_SIM22 = imwarp(SIM22,affine2d(tform2'),'OutputView',Rfixed);

%% rot270
SIM10 = rot90(SIM10,3);
SIM11 = rot90(SIM11,3);
SIM12 = rot90(SIM12,3);

r_SIM00 = rot90(r_SIM00,3);
r_SIM01 = rot90(r_SIM01,3);
r_SIM02 = rot90(r_SIM02,3);

r_SIM20 = rot90(r_SIM20,3);
r_SIM21 = rot90(r_SIM21,3);
r_SIM22 = rot90(r_SIM22,3);

imgch0= r_SIM00+r_SIM01+r_SIM02;
imgch1= SIM10+SIM11+SIM12;
imgch2= r_SIM20+r_SIM21+r_SIM22; 

%% output


name_SIM10= fullfile(folder, strcat('reg2_SIM10-',num2str(n),'.tif'))
SupFun.tiffWrite(SIM10,name_SIM10)

name_SIM11= fullfile(folder, strcat('reg2_SIM11-',num2str(n),'.tif'))
SupFun.tiffWrite(SIM11,name_SIM11)

name_SIM12= fullfile(folder, strcat('reg2_SIM12-',num2str(n),'.tif'))
SupFun.tiffWrite(SIM12,name_SIM12)

name_SIM00= fullfile(folder, strcat('reg2_SIM00-',num2str(n),'.tif'));
SupFun.tiffWrite(r_SIM00,name_SIM00)

name_SIM01= fullfile(folder, strcat('reg2_SIM01-',num2str(n),'.tif'));
SupFun.tiffWrite(r_SIM01,name_SIM01)

name_SIM02= fullfile(folder, strcat('reg2_SIM02-',num2str(n),'.tif'));
SupFun.tiffWrite(r_SIM02,name_SIM02)

name_SIM20= fullfile(folder, strcat('reg2_SIM20-',num2str(n),'.tif'));
SupFun.tiffWrite(r_SIM20,name_SIM20)

name_SIM21= fullfile(folder, strcat('reg2_SIM21-',num2str(n),'.tif'));
SupFun.tiffWrite(r_SIM21,name_SIM21)

name_SIM22= fullfile(folder, strcat('reg2_SIM22-',num2str(n),'.tif'));
SupFun.tiffWrite(r_SIM22,name_SIM22)
    

toc
end
close all
toc