

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               Slice by slice SIM processing
%
%               -edge mirroring added
%               -Unmixing at native resolution, then upsampling via zero
%               padding
%               -throwing out slices that have bad condition number of
%               unmixing matrix
%
%
%               Reto Fiolka 15th September 2021
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc


%% file paths:

DatFolder='/archive/bioinformatics/Danuser_lab/Fiolka/Manuscripts/2021_OPMSIM/MatlabCode/ForNatureMethods/OPSIM_ExampleData/'



%% general settings
dpix=6.5/(57);     % pixel size in microns
per=0.37%*561/488;  % expected period for SIM pattern, for 488nm exc, about 0.370microns
Res=0.280          % expected OPM raw resolution, defines cut off frequency (1/Res).
mod=5.2 %estimated modulation depth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters and options

saving=1 % save data when 1 is selected. Select 0 to plot intermediate data


       
multiC=0  %option for multicolor processing: for second stack, expected period will get reduced by factor 488/561

preRL=1   %RL Deconvolution before SIM processing
Zoff=0    %z shift of direction 1
bd=10     % border around raw images

dampen=1  % use edge mirroring to reduce gibbs' oscillations
w=20;     % width for edge mirroring

lu=[10:180];%[65:105];  %Selected the slices to be processed. For a range use e.g. lu=[50:70]

t=[0];

susli=0;    % sum adjacent slices before processing


miCo=1  %% missing cone infill when 1 is selected
CenterWeight=0.15;  %factor to multiply centeral missing cone in-fill.

NotchW=1/40;      %Notch filter width, as fraction of cut-off freqeuncy

Nweight=0.3;      % weight of notch in-fill


apod3D=0
%conT=1;   %threshold for condition number; if it is exceeded, use MIP estimate for phase and peak position



if saving==1
    plo=0
else
    plo=1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create directories to save dat%PSF_folder='/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPMSIM/BeadsSIM/CoverSlipBeadsPSF'a

if saving
    recFolder=[DatFolder '/SIM2/']%PSF_folder='/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPMSIM/BeadsSIM/CoverSlipBeadsPSF'
   mkdir(recFolder)
   recFolderS=[DatFolder '/SIM2/Stacks']
   mkdir(recFolderS)%PSF_folder='/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPMSIM/BeadsSIM/CoverSlipBeadsPSF'
   recFolderM=[DatFolder '/SIM2/MIPS']
   mkdir(recFolderM)
end%PSF_folder='/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPMSIM/BeadsSIM/CoverSlipBeadsPSF'


%% read PSF

PSF=double(imread('/archive/bioinformatics/Danuser_lab/Fiolka/Manuscripts/2021_OPMSIM/MatlabCode/ForNatureMethods/OPSIM_ExampleData/PSF_NA1p1_520nm_x100nm.tif'));
PSFsim=double(imread('/archive/bioinformatics/Danuser_lab/Fiolka/Manuscripts/2021_OPMSIM/MatlabCode/ForNatureMethods/OPSIM_ExampleData/SIMpsf.tif'));
PSF=PSF(31:end-31,31:end-31);

for n=t
       if n==1 && multiC==1   %option for multicolor processing
           per=0.38*561/488;
       end
% 
SIM00=double(SupFun.stackRead([DatFolder 'reg2_SIM00-' num2str(n) '.tif']));
SIM01=double(SupFun.stackRead([DatFolder 'reg2_SIM01-' num2str(n) '.tif']));
SIM02=double(SupFun.stackRead([DatFolder 'reg2_SIM02-' num2str(n) '.tif']));


SIM10=double(SupFun.stackRead([DatFolder 'reg2_SIM10-' num2str(n) '.tif']));
SIM11=double(SupFun.stackRead([DatFolder 'reg2_SIM11-' num2str(n) '.tif']));
SIM12=double(SupFun.stackRead([DatFolder 'reg2_SIM12-' num2str(n) '.tif']));

SIM20=double(SupFun.stackRead([DatFolder 'reg2_SIM20-' num2str(n) '.tif']));
SIM21=double(SupFun.stackRead([DatFolder 'reg2_SIM21-' num2str(n) '.tif']));
SIM22=double(SupFun.stackRead([DatFolder 'reg2_SIM22-' num2str(n) '.tif']));


[ydim,xdim,zdim]=size(SIM00);

%% create a square x-y ROI, using the smaller of the two lateral dimensions
if xdim<ydim      
    xdim=xdim;
elseif ydim<xdim
    xdim=ydim;
else
    xdim=xdim;
    
end



SIM00=SIM00(bd:xdim-bd,bd:xdim-bd,:);
SIM01=SIM01(bd:xdim-bd,bd:xdim-bd,:);
SIM02=SIM02(bd:xdim-bd,bd:xdim-bd,:);

SIM10=SIM10(bd:xdim-bd,bd:xdim-bd,:);
SIM11=SIM11(bd:xdim-bd,bd:xdim-bd,:);
SIM12=SIM12(bd:xdim-bd,bd:xdim-bd,:);

%%%%%%%%%%

    SIM10=circshift(SIM10,Zoff,3);
       SIM11=circshift(SIM11,Zoff,3);
           SIM12=circshift(SIM12,Zoff,3);
%%%%%%%%%%%%%%%

SIM20=SIM20(bd:xdim-bd,bd:xdim-bd,:);
SIM21=SIM21(bd:xdim-bd,bd:xdim-bd,:);
SIM22=SIM22(bd:xdim-bd,bd:xdim-bd,:);



[ydim,xdim,zdim]=size(SIM00);

%%  Compute pattern and cutoff frequency based on xy pixel size

Kc=1/(Res/dpix)*xdim;     % cut off frequency based on estimated OPM resolution

%% compute annulus around expected pattern frequency

freq=1/(1.07*per/dpix)*xdim;
block=SupFun.kreis(freq,xdim,xdim,0,0);
block2=SupFun.kreis(freq*1.11,xdim,xdim,0,0);
iblock=find(block(:)==1);

test=block2-block;
iblock2=find(test(:)==0);
fla0=0;fla1=0;fla2=0;


%% estimate pattern and phase from MIPs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first direction, Dir0:

miSIM00=max(SIM00,[],3);miSIM01=max(SIM01,[],3);miSIM02=max(SIM02,[],3);

bl=size(miSIM00,1);bl2=round(bl/2);
  
    
    AFFT=abs(fftshift(fft2(miSIM00)))+abs(fftshift(fft2(miSIM01)))+abs(fftshift(fft2(miSIM02)));
    AFFT(iblock2)=0;
    AFFT(bl2-10:bl2+10,:)=0;
    AFFT(1:bl2,bl2:end)=0;
    AFFT(bl2:end,1:bl2)=0;
    AFFT(:,bl2-10:bl2+10)=0;
    [afft,b]=size(AFFT);
    if plo
    figure;imshow(abs(AFFT),[])
    end
    indmax= find(imregionalmax(AFFT) == 1);        %indexes of maxima
    [peaks,sortind]= sort(AFFT(indmax),1,'descend');
    i1=1; i2=2;
    [ky1,kx1]=ind2sub(size(AFFT),indmax(sortind(i1)));
    [ky2,kx2]=ind2sub(size(AFFT),indmax(sortind(i2)));
  
    kx11=kx1;ky11=ky1;
    kx22=kx2;ky22=ky2;
    
    u1=max([kx11,ky11]);
    u2=max([kx22,ky22]);
    if u1>u2;
        dkx=(kx11-kx22)/2;
        dky=(ky11-ky22)/2;
    elseif u1<u2
        dkx=(kx22-kx11)/2;

        dky=(ky22-ky11)/2;
    elseif u1==u2
        dkx=(max([kx11,kx22])-min([kx11,kx22]))/2;
        dky=(max([ky11,ky22])-min([ky11,ky22]))/2;
    end
    dkx0=dkx;dky0=dky;
AFFT0=(fftshift(fft2(miSIM00)));
AFFT1=(fftshift(fft2(miSIM01)));
AFFT2=(fftshift(fft2(miSIM02)));

phasex01=angle(conj(AFFT0(ky2,kx2))*AFFT1(ky2,kx2));
phasex02=angle(conj(AFFT1(ky2,kx2))*AFFT2(ky2,kx2));
pha1=phasex01*180/pi
pha2=phasex02*180/pi

%% DIR1
miSIM10=max(SIM10,[],3);miSIM11=max(SIM11,[],3);miSIM12=max(SIM12,[],3);


    AFFT=abs(fftshift(fft2(miSIM10)))+abs(fftshift(fft2(miSIM11)))+abs(fftshift(fft2(miSIM12)));
    AFFT(iblock2)=0;
    AFFT(bl2-10:bl2+10,:)=0;
    [afft,b]=size(AFFT);
    if plo
    figure;imshow(abs(AFFT),[])
    end
    indmax= find(imregionalmax(AFFT) == 1);        %indexes of maxima
    [peaks,sortind]= sort(AFFT(indmax),1,'descend');
    i1=1; i2=2;
    [ky1,kx1]=ind2sub(size(AFFT),indmax(sortind(i1)));
    [ky2,kx2]=ind2sub(size(AFFT),indmax(sortind(i2)));
  
    kx11=kx1;ky11=ky1;
    kx22=kx2;ky22=ky2;
    
    u1=max([kx11,ky11]);
    u2=max([kx22,ky22]);
    if u1>u2;
        dkx=(kx11-kx22)/2;
        dky=(ky11-ky22)/2;
    elseif u1<u2
        dkx=(kx22-kx11)/2;

        dky=(ky22-ky11)/2;
    elseif u1==u2
        dkx=(max([kx11,kx22])-min([kx11,kx22]))/2;
        dky=(max([ky11,ky22])-min([ky11,ky22]))/2;
    end
    dkx1=dkx;dky1=dky;
AFFT0=(fftshift(fft2(miSIM10)));
AFFT1=(fftshift(fft2(miSIM11)));
AFFT2=(fftshift(fft2(miSIM12)));

phasex11=angle(conj(AFFT0(ky2,kx2))*AFFT1(ky2,kx2));
phasex12=angle(conj(AFFT1(ky2,kx2))*AFFT2(ky2,kx2));
pha1=phasex11*180/pi
pha2=phasex12*180/pi

%% DIR 2
miSIM20=max(SIM20,[],3);miSIM21=max(SIM21,[],3);miSIM22=max(SIM22,[],3);


    AFFT=abs(fftshift(fft2(miSIM20)))+abs(fftshift(fft2(miSIM21)))+abs(fftshift(fft2(miSIM22)));
     AFFT(iblock2)=0;
      AFFT(bl2-10:bl2+10,:)=0;
    AFFT(:,bl2-10:bl2+10)=0;
    
       AFFT(1:bl2,1:bl2)=0;
      AFFT(bl2:end,bl2:end)=0;
      [afft,b]=size(AFFT);
    [afft,b]=size(AFFT);
    if plo
    figure;imshow(abs(AFFT),[])
    end
    indmax= find(imregionalmax(AFFT) == 1);        %indexes of maxima
    [peaks,sortind]= sort(AFFT(indmax),1,'descend');
    i1=1; i2=2;
    [ky1,kx1]=ind2sub(size(AFFT),indmax(sortind(i1)));
    [ky2,kx2]=ind2sub(size(AFFT),indmax(sortind(i2)));
  
    kx11=kx1;ky11=ky1;
    kx22=kx2;ky22=ky2;
    
    u1=max([kx11,ky11]);
    u2=max([kx22,ky22]);
    if u1>u2;
        dkx=(kx11-kx22)/2;
        dky=(ky11-ky22)/2;
    elseif u1<u2
        dkx=(kx22-kx11)/2;

        dky=(ky22-ky11)/2;
    elseif u1==u2
        dkx=(max([kx11,kx22])-min([kx11,kx22]))/2;
        dky=(max([ky11,ky22])-min([ky11,ky22]))/2;
    end
    dkx2=dkx;dky2=dky;
AFFT0=(fftshift(fft2(miSIM20)));
AFFT1=(fftshift(fft2(miSIM21)));
AFFT2=(fftshift(fft2(miSIM22)));

phasex21=angle(conj(AFFT0(ky2,kx2))*AFFT1(ky2,kx2));
phasex22=angle(conj(AFFT1(ky2,kx2))*AFFT2(ky2,kx2));
pha1=phasex22*180/pi
pha2=phasex22*180/pi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% start processing one slice by slice

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count=1; % book keeping index starting at 1

for l=lu
  
    

X00=(squeeze(SIM00(:,:,l)));
X01=(squeeze(SIM01(:,:,l)));
X02=(squeeze(SIM02(:,:,l)));

X10=(squeeze(SIM10(:,:,l)));
X11=(squeeze(SIM11(:,:,l)));
X12=(squeeze(SIM12(:,:,l)));

X20=(squeeze(SIM20(:,:,l)));
X21=(squeeze(SIM21(:,:,l)));
X22=(squeeze(SIM22(:,:,l)));



%% edge dampening option
if dampen

X00=SupFun.DampEdgeE(X00,w);
X01=SupFun.DampEdgeE(X01,w);
X02=SupFun.DampEdgeE(X02,w);

X10=SupFun.DampEdgeE(X10,w);
X11=SupFun.DampEdgeE(X11,w);
X12=SupFun.DampEdgeE(X12,w);

X20=SupFun.DampEdgeE(X20,w);
X21=SupFun.DampEdgeE(X21,w);
X22=SupFun.DampEdgeE(X22,w);
end

%% SIM processing of Dir 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%% estimate pattern frequency from peak   
    temp=X00;
    AFFT=abs(fftshift(fft2(temp)));
    AFFT(iblock2)=0;
    AFFT(bl2-10:bl2+10,:)=0;
      AFFT(1:bl2,bl2:end)=0;
      AFFT(bl2:end,1:bl2)=0;
    AFFT(:,bl2-10:bl2+10)=0;
    [afft,b]=size(AFFT);
    if plo
    figure;imshow(abs(AFFT),[])
    end
    indmax= find(imregionalmax(AFFT) == 1);        %indexes of maxima
    [peaks,sortind]= sort(AFFT(indmax),1,'descend');
 
    i1=1; i2=2;
    [ky1,kx1]=ind2sub(size(AFFT),indmax(sortind(i1)));
    [ky2,kx2]=ind2sub(size(AFFT),indmax(sortind(i2)));

    
    kx11=kx1;ky11=ky1;
    kx22=kx2;ky22=ky2;
    
    u1=max([kx11,ky11]);
    u2=max([kx22,ky22]);
    if u1>u2;
        dkx=(kx11-kx22)/2;
        dky=(ky11-ky22)/2;
    elseif u1<u2
        dkx=(kx22-kx11)/2;

        dky=(ky22-ky11)/2;
    elseif u1==u2
        dkx=(max([kx11,kx22])-min([kx11,kx22]))/2;
        dky=(max([ky11,ky22])-min([ky11,ky22]))/2;
    end
    
    Tx=afft*dpix/dkx;
    Ty=afft*dpix/dky;
    K=sqrt(dkx^2+dky^2);
    T=afft*dpix/K 


%% determine phase shifts

[xdim,ydim]=size(X00);

AFFT0=(fftshift(fft2(X00)));
AFFT1=(fftshift(fft2(X01)));
AFFT2=(fftshift(fft2(X02)));

phasex1=angle(conj(AFFT0(ky2,kx2))*AFFT1(ky2,kx2));
phasex2=angle(conj(AFFT1(ky2,kx2))*AFFT2(ky2,kx2));
pha1=phasex1*180/pi
pha2=phasex2*180/pi

%% RL deconvolution of raw data before processing
if preRL
X00=deconvlucy(X00,PSF,10);
X01=deconvlucy(X01,PSF,10);
X02=deconvlucy(X02,PSF,10);
end
bright0=X00+X01+X02;        %% deconvolved widefield image by summing images from direction 0

XF00=fftshift(fft2(X00));
XF01=fftshift(fft2(X01));
XF02=fftshift(fft2(X02));


%% unmixing of bands using estimated phase steps
% phasex1=120/180*pi;phasex2=phasex1;
%   [XF00,XF01,XF02,M]=Demod(XF00,XF01,XF02,0,phasex1,phasex2+phasex1);
%  co=cond(M) 
%  if co>conT; 

[XF00,XF01,XF02,M]=SupFun.Demod(XF00,XF01,XF02,0,phasex01,phasex01+phasex02);
% co=cond(M)
%  disp('using MIP phase values')
%  end

disp('Demodulation of spectra finished')

%% zero pad spectra 
XF00=padarray(XF00,[round(size(XF00,2)/2),round(size(XF00,2)/2)]);
XF01=padarray(XF01,[round(size(XF01,2)/2),round(size(XF01,2)/2)]);
XF02=padarray(XF02,[round(size(XF02,2)/2),round(size(XF02,2)/2)]);

[xdim,ydim]=size(XF00);


%% shifting of the spectra in real space
% needs work, circshift is only integer pixels, but somehow works
% acceptably


S00f=XF00;

S01f=circshift(XF01,[-round(dky0),-round(dkx0)]);
S02f=circshift(XF02,[round(dky0),round(dkx0)]);  
disp('shift using MIP values')


%% Determination of the initial phase
d2x=round(xdim/2);
d2y=d2x
d=60;
phas0=S00f(d2y-d:d2y+d,d2x-d:d2x+d);phas0(d2y-3:d2y+3,d2x-3:d2x+3)=0;
phas1=S01f(d2y-d:d2y+d,d2x-d:d2x+d);phas1(d2y-3:d2y+3,d2x-3:d2x+3)=0;%    
phas2=S02f(d2y-d:d2y+d,d2x-d:d2x+d);phas2(d2y-3:d2y+3,d2x-3:d2x+3)=0;

clear i;
clear m;

%difference in the phase angles:

c=sum(sum(conj(phas0).*phas1));
phaseave1=angle(c);

c=sum(sum(conj(phas0).*phas2));
phaseave2=angle(c);
mean(phaseave1*180/pi);
mean(phaseave2*180/pi);

S01f=S01f*exp(-sqrt(-1)*phaseave1);  % correction of global phase angle
S02f=S02f*exp(-sqrt(-1)*phaseave2);


%% Apodization masks, missing cone in-fill and notch filter

mC00=SupFun.kreis(round(Kc/4),xdim,xdim,0,0);             % Missing cone infill
mC01=SupFun.kreis(round(Kc*NotchW),xdim,xdim,-dkx0,-dky0);  % Notch filter mask 1
mC02=SupFun.kreis(round(Kc*NotchW),xdim,xdim,dkx0,dky0);    % Notch filter mask 2

Apo00=SupFun.kreis(round(Kc),xdim,xdim,0,0);       % binary apodization masks
Apo01=SupFun.kreis(round(Kc),xdim,xdim,-dkx0,-dky0); % binary apodization masks
Apo02=SupFun.kreis(round(Kc),xdim,xdim,dkx0,dky0);   % binary apodization masks

S00f=S00f.*Apo00;S01f=S01f.*Apo01;S02f=S02f.*Apo02;  %binary apodization of bands

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SIM processing for Dir1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% estimate pattern frequency    
    temp=X11;
    AFFT=abs(fftshift(fft2(temp)));
    AFFT(iblock2)=0;
    AFFT(bl2-10:bl2+10,:)=0;
    [afft,b]=size(AFFT);
    if plo
    figure;imshow(abs(AFFT),[])
    end
    indmax= find(imregionalmax(AFFT) == 1);        %indexes of maxima
    [peaks,sortind]= sort(AFFT(indmax),1,'descend');
    i1=1; i2=2;
    [ky1,kx1]=ind2sub(size(AFFT),indmax(sortind(i1)));
    [ky2,kx2]=ind2sub(size(AFFT),indmax(sortind(i2)));

    
    kx11=kx1;ky11=ky1;
    kx22=kx2;ky22=ky2;
    
    u1=max([kx11,ky11]);
    u2=max([kx22,ky22]);
    if u1>u2;
        dkx=(kx11-kx22)/2;
        dky=(ky11-ky22)/2;
    elseif u1<u2
        dkx=(kx22-kx11)/2;

        dky=(ky22-ky11)/2;
    elseif u1==u2
        dkx=(max([kx11,kx22])-min([kx11,kx22]))/2;
        dky=(max([ky11,ky22])-min([ky11,ky22]))/2;
    end
    
    Tx=afft*dpix/dkx;
    Ty=afft*dpix/dky;
    K=sqrt(dkx^2+dky^2);
    T=afft*dpix/K


%% determine phase shifts:

AFFT0=(fftshift(fft2(X10)));
AFFT1=(fftshift(fft2(X11)));
AFFT2=(fftshift(fft2(X12)));

phasex1=angle(conj(AFFT0(ky2,kx2))*AFFT1(ky2,kx2));
phasex2=angle(conj(AFFT1(ky2,kx2))*AFFT2(ky2,kx2));
pha1=phasex1*180/pi;
pha2=phasex2*180/pi;

if preRL
X10=deconvlucy(X10,PSF,10);
X11=deconvlucy(X11,PSF,10);
X12=deconvlucy(X12,PSF,10);
end

bright1=X00+X01+X02;
XF10=fftshift(fft2(X10));
XF11=fftshift(fft2(X11));
XF12=fftshift(fft2(X12));

%  phasex1=120*pi/180;
%  phasex2=phasex1;
%  [XF10,XF11,XF12,M]=Demod(XF10,XF11,XF12,0,phasex1,phasex2+phasex1);
%  co=cond(M); if co>conT; 
[XF10,XF11,XF12,M]=SupFun.Demod(XF10,XF11,XF12,0,phasex11,phasex12+phasex11);
%  disp('using MIP determined values')
% end
disp('Demodulation of spectra finished')

%% zero pad spectra
XF10=padarray(XF10,[round(size(XF10,2)/2),round(size(XF10,2)/2)]);
XF11=padarray(XF11,[round(size(XF11,2)/2),round(size(XF11,2)/2)]);
XF12=padarray(XF12,[round(size(XF12,2)/2),round(size(XF12,2)/2)]);
[xdim,ydim]=size(XF00);

%% shifting of the spectra in real space

[px,py] = meshgrid(0:afft-1);

S10f=XF10;


S11f=circshift(XF11,[-round(dky1),-round(dkx1)]);
S12f=circshift(XF12,[round(dky1),round(dkx1)]);
disp('using shift values from MIP')

 

%% Determination of the initial phase

d=60;
phas0=S10f(d2y-d:d2y+d,d2x-d:d2x+d);phas0(d2y-3:d2y+3,d2x-3:d2x+3)=0;
phas1=S11f(d2y-d:d2y+d,d2x-d:d2x+d);phas1(d2y-3:d2y+3,d2x-3:d2x+3)=0;%
phas2=S12f(d2y-d:d2y+d,d2x-d:d2x+d);phas2(d2y-3:d2y+3,d2x-3:d2x+3)=0;

clear i;
clear m;

%difference in the phase angles:

c=sum(sum(conj(phas0).*phas1));
phaseave1=angle(c);

c=sum(sum(conj(phas0).*phas2));
phaseave2=angle(c);
mean(phaseave1*180/pi);
mean(phaseave2*180/pi);

m1=sum(sum(abs(phas0))/sum(abs(phas1))); 
m2=sum(sum(abs(phas0))/sum(abs(phas2)));

S11f=S11f*exp(-sqrt(-1)*phaseave1);
S12f=S12f*exp(-sqrt(-1)*phaseave2);


%% Apodization

mC10=SupFun.kreis(round(Kc/4),xdim,xdim,0,0);
mC11=SupFun.kreis(round(Kc*NotchW),xdim,xdim,-dkx1,-dky1);
mC12=SupFun.kreis(round(Kc*NotchW),xdim,xdim,dkx1,dky1);

Apo10=SupFun.kreis(round(Kc),xdim,xdim,0,0);
Apo11=SupFun.kreis(round(Kc),xdim,xdim,-dkx1,-dky1);
Apo12=SupFun.kreis(round(Kc),xdim,xdim,dkx1,dky1);


S10f=S10f.*Apo10;S11f=S11f.*Apo11;S12f=S12f.*Apo12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SIM processing of Dir2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% estimate pattern frequency    
    temp=X20;
 
    AFFT=abs(fftshift(fft2(temp)));
    AFFT(iblock2)=0;
    AFFT(bl2-10:bl2+10,:)=0;
    AFFT(:,bl2-10:bl2+10)=0;
    
      AFFT(1:bl2,1:bl2)=0;
      AFFT(bl2:end,bl2:end)=0;
      [afft,b]=size(AFFT);
      if plo
      figure;imshow(abs(AFFT),[])
      end
    indmax= find(imregionalmax(AFFT) == 1);        %indexes of maxima
    [peaks,sortind]= sort(AFFT(indmax),1,'descend');

    i1=1; i2=2;
    [ky1,kx1]=ind2sub(size(AFFT),indmax(sortind(i1)));
    [ky2,kx2]=ind2sub(size(AFFT),indmax(sortind(i2)));

    
    kx11=kx1;ky11=ky1;
    kx22=kx2;ky22=ky2;
    
    u1=max([kx11,ky11]);
    u2=max([kx22,ky22]);
    if u1>u2;
        dkx=(kx11-kx22)/2;
        dky=(ky11-ky22)/2;
    elseif u1<u2
        dkx=(kx22-kx11)/2;

        dky=(ky22-ky11)/2;
    elseif u1==u2
        dkx=(max([kx11,kx22])-min([kx11,kx22]))/2;
        dky=(max([ky11,ky22])-min([ky11,ky22]))/2;
    end
    
    Tx=afft*dpix/dkx;
    Ty=afft*dpix/dky;
    K=sqrt(dkx^2+dky^2);
    T=afft*dpix/K


%% determine phase shifts:

AFFT0=(fftshift(fft2(X20)));
AFFT1=(fftshift(fft2(X21)));
AFFT2=(fftshift(fft2(X22)));

phasex1=angle(conj(AFFT0(ky2,kx2))*AFFT1(ky2,kx2));
phasex2=angle(conj(AFFT1(ky2,kx2))*AFFT2(ky2,kx2));
pha1=phasex1*180/pi;
pha2=phasex2*180/pi;

if preRL
X20=deconvlucy(X20,PSF,10);
X21=deconvlucy(X21,PSF,10);
X22=deconvlucy(X22,PSF,10);
end

bright2=X20+X21+X22;
XF20=fftshift(fft2(X20));
XF21=fftshift(fft2(X21));
XF22=fftshift(fft2(X22));

%% unmix bands
%  phasex1=120*pi/180;
%  phasex2=phasex1;
%  [XF20,XF21,XF22,M]=Demod(XF20,XF21,XF22,0,phasex1,phasex2+phasex1);
% co=cond(M); if co>conT; 
[XF20,XF21,XF22,M]=SupFun.Demod(XF20,XF21,XF22,0,phasex21,phasex22+phasex21);
% disp('using MIP determined values')
%  end

disp('Demodulation of spectra finished')

%% zero pad spectra
XF20=padarray(XF20,[round(size(XF20,2)/2),round(size(XF20,2)/2)]);
XF21=padarray(XF21,[round(size(XF21,2)/2),round(size(XF21,2)/2)]);
XF22=padarray(XF22,[round(size(XF22,2)/2),round(size(XF22,2)/2)]);
[xdim,ydim]=size(XF00);

%% shifting of the spectra in real space

S20f=XF20;

S21f=circshift(XF21,[-round(dky2),-round(dkx2)]);
S22f=circshift(XF22,[round(dky2),round(dkx2)]);
disp('shift using MIP values')

 
%% Determination of the initial phase

d=60;
phas0=S20f(d2y-d:d2y+d,d2x-d:d2x+d);phas0(d2y-3:d2y+3,d2x-3:d2x+3)=0;
phas1=S21f(d2y-d:d2y+d,d2x-d:d2x+d);phas1(d2y-3:d2y+3,d2x-3:d2x+3)=0;
phas2=S22f(d2y-d:d2y+d,d2x-d:d2x+d);phas2(d2y-3:d2y+3,d2x-3:d2x+3)=0;

clear i;
clear m;

%difference in the phase angles:

c=sum(sum(conj(phas0).*phas1));
phaseave1=angle(c);

c=sum(sum(conj(phas0).*phas2));
phaseave2=angle(c);
mean(phaseave1*180/pi);
mean(phaseave2*180/pi);

m1=sum(sum(abs(phas0))/sum(abs(phas1))); 
m2=sum(sum(abs(phas0))/sum(abs(phas2)));

S21f=S21f*exp(-sqrt(-1)*phaseave1);
S22f=S22f*exp(-sqrt(-1)*phaseave2);


%% Apodization

mC20=SupFun.kreis(round(Kc/4),xdim,xdim,0,0);
mC21=SupFun.kreis(round(Kc*NotchW),xdim,xdim,-dkx2,-dky2);
mC22=SupFun.kreis(round(Kc*NotchW),xdim,xdim,dkx2,dky2);

Apo20=SupFun.kreis(round(Kc),xdim,xdim,0,0);
Apo21=SupFun.kreis(round(Kc),xdim,xdim,-dkx2,-dky2);
Apo22=SupFun.kreis(round(Kc),xdim,xdim,dkx2,dky2);


S20f=S20f.*Apo20;S21f=S21f.*Apo21;S22f=S22f.*Apo22;



%% Assemble spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% compute weights of different bands
%disp('mod value')
%mod
ApoTot=Apo00+mod*Apo01+mod*Apo02+Apo10+mod*Apo11+mod*Apo12+Apo20+mod*Apo21+mod*Apo22;
indS=find(ApoTot==5.2);

ApoTot(indS)=mod*3.5;

inda=find(ApoTot>0);

% sum all bands together
Recon=(S01f+S02f+S00f+S11f+S12f+S10f+S21f+S22f+S20f);
Recon(inda)=Recon(inda)./ApoTot(inda);


% find indices of missing cone in-fill
ind00=find(mC00==1);
ind01=find(mC01==1);
ind02=find(mC02==1);

ind10=find(mC10==1);
ind11=find(mC11==1);
ind12=find(mC12==1);

ind20=find(mC20==1);
ind21=find(mC21==1); 
ind22=find(mC22==1);

% replace missing cone region with sideband information
if miCo
Recon(ind00)=(S01f(ind00)+S11f(ind00)+S21f(ind00)+S22f(ind00)+S12f(ind00)+S02f(ind00))*CenterWeight;
end

% notch filter / repace pattern peaks with zero order band information
Recon(ind01)=S00f(ind01)*Nweight;
Recon(ind02)=S00f(ind02)*Nweight;
Recon(ind11)=S00f(ind11)*Nweight;
Recon(ind12)=S00f(ind12)*Nweight;
Recon(ind21)=S00f(ind21)*Nweight;
Recon(ind22)=S00f(ind22)*Nweight;
%% final RL deconvolution of SIM reconstruction


    
Recon2=real(ifft2(ifftshift(Recon)));
Recon3=deconvlucy(Recon2,PSFsim,2);
if plo
RLRecon=(fftshift(fft2(Recon3)));
figure; imshow(log(abs(Recon)),[])
figure; imshow(log(abs(RLRecon)),[])
figure;imshow(Recon3,[])
end



RecStack(:,:,count)=Recon3;
WFSum(:,:,count)=bright0+bright1+bright2;  %summed WF image, deconvolved




disp(['plane ' num2str(count) ' processed'])
count=count+1

end

if apod3D
    
    T=SupFun.Triangle3D(1*round(K+Kc),xdim,xdim,size(RecStack,3),0,0,0);
    
        fusion=(fftshift(fftn(RecStack)));
    
    fusion=fusion.*T;
    RecStack=real(ifftn(ifftshift(fusion)));
    disp('done with 3D apodization')
    
end


%% saving of data
if saving

%normalize stack    

RecStack=RecStack/max(RecStack(:));

WFSum=WFSum/max(WFSum(:));

RecStack=RecStack*20000;

WFSum=WFSum*20000;


SupFun.stackWrite(uint16(WFSum),[recFolder 'Widefieldstack-' num2str(n) '.tif']);

SupFun.stackWrite(uint16(RecStack),[recFolderS '/SIMstackPlus-' num2str(n) '.tif']);
temp=max(RecStack,[],3); imwrite(uint16(temp), [recFolderM '/SIMstackPlusMIP-' num2str(n) '.tif']);

end


end

