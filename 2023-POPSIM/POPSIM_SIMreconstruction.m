

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               Projective Oblique Structured illumination Microscopy processing
%
%               -edge mirroring added
%               -Unmixing at native resolution, then upsampling via zero
%               padding
%               -throwing out slices that have bad condition number of
%               unmixing matrix
%
%
%               Reto Fiolka 23rd August 2023
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc


%% file path:

DatFolder='/archive/bioinformatics/Danuser_lab/Fiolka/Manuscripts/POPSIM/CodeToShare/ExampleData/'

%% settings
dpix=6.5/(57);     % pixel size in microns
per=0.365 % expected period for SIM pattern, for 488nm exc, about 0.370microns
Res=0.290 % expected OPM raw resolution, defines cut off frequency (1/Res).


saving=0 % save data if saving=1

t=[1]; % selected timepoint
   
preRL=1   %RL Deconvolution before SIM processing

RL1=10; %number of iterations for RL deconvolution.

miCo=1  %% missing cone infill
MiCoRa=1/3; %Radius of missing cone in-fill, in relation to cut-off frequency
CenterWeight=0.07;  %factor to multiply centeral missing cone in-fill.
WFweight=0.05          % how much DC band is mixed into central region

NotchW=1/50;      %Notch filter width, as fraction of cut-off freqeuncy
Nweight=0.0;      % weight of notch in-fill

apod=1%% triangular apodization option
scaleT=1; %scale of the triangular cone


dampen=1  % use edge mirroring to reduce gibbs' oscillations
w=10;     % width for edge mirroring
db=5     % border of image that will be truncated

if saving==1
    plo=0;
else
    plo=1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create directories to save data

if saving
   recFolder=[DatFolder '/SIM/']
   mkdir(recFolder)

   recFolderM=[DatFolder '/SIM/Reconstuctions']
   mkdir(recFolderM)
end



%% read PSF

PSF=double(imread([DatFolder '/PSF_NA1p1_520nm_x100nm.tif']));
PSF=PSF(31:end-31,31:end-31);
PSFsim=double(imread([DatFolder '/SIMpsf.tif']));



for n=t(1):t(end)


SIM00=double(imread([DatFolder 'reg2_SIM00-' num2str(n) '.tif']));
SIM01=double(imread([DatFolder 'reg2_SIM01-' num2str(n) '.tif']));
SIM02=double(imread([DatFolder 'reg2_SIM02-' num2str(n) '.tif']));


SIM10=double(imread([DatFolder 'reg2_SIM10-' num2str(n) '.tif']));
SIM11=double(imread([DatFolder 'reg2_SIM11-' num2str(n) '.tif']));
SIM12=double(imread([DatFolder 'reg2_SIM12-' num2str(n) '.tif']));

SIM20=double(imread([DatFolder 'reg2_SIM20-' num2str(n) '.tif']));
SIM21=double(imread([DatFolder 'reg2_SIM21-' num2str(n) '.tif']));
SIM22=double(imread([DatFolder 'reg2_SIM22-' num2str(n) '.tif']));

WFSum=SIM10+SIM11+SIM12+SIM00+SIM01+SIM02+SIM20+SIM21+SIM22;


[ydim,xdim]=size(SIM00);

%% create a square x-y ROI, using the smaller of the two lateral dimensions

if xdim<ydim      
    xdim=xdim;
elseif ydim<xdim
    xdim=ydim;
else
    xdim=xdim;
    
end

SIM00=SIM00(db:xdim-db,db:xdim-db);
SIM01=SIM01(db:xdim-db,db:xdim-db);
SIM02=SIM02(db:xdim-db,db:xdim-db);

SIM10=SIM10(db:xdim-db,db:xdim-db);
SIM11=SIM11(db:xdim-db,db:xdim-db);
SIM12=SIM12(db:xdim-db,db:xdim-db);


SIM20=SIM20(db:xdim-db,db:xdim-db);
SIM21=SIM21(db:xdim-db,db:xdim-db);
SIM22=SIM22(db:xdim-db,db:xdim-db);

WFSum=WFSum(db:xdim-db,db:xdim-db);

[ydim,xdim]=size(SIM00);


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


%% estimate pattern and phase 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Dir0:

miSIM00=(SIM00);miSIM01=SIM01;miSIM02=SIM02;

bl=size(miSIM00,1);bl2=round(bl/2);
    
    AFFT=abs(fftshift(fft2(miSIM00)))+abs(fftshift(fft2(miSIM01)))+abs(fftshift(fft2(miSIM02)));
 
test=SupFun.kreis(6,size(AFFT,1),size(AFFT,2),sind(60)*freq,-cosd(60)*freq)+SupFun.kreis(6,size(AFFT,1),size(AFFT,2),-sind(60)*freq,cosd(60)*freq);


AFFT=AFFT.*test;
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
    Tx=afft*dpix/dkx;
    Ty=afft*dpix/dky;
    K=sqrt(dkx^2+dky^2);
    T=afft*dpix/K

%% DIR1
miSIM10=SIM10;miSIM11=SIM11;miSIM12=SIM12;

    AFFT=abs(fftshift(fft2(miSIM10)))+abs(fftshift(fft2(miSIM11)))+abs(fftshift(fft2(miSIM12)));
    
    test=SupFun.kreis(10,size(AFFT,1),size(AFFT,2),0,-freq*1.02)+SupFun.kreis(10,size(AFFT,1),size(AFFT,2),0,freq*1.02);
  
    AFFT=AFFT.*test;
 
    
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
    Tx=afft*dpix/dkx;
    Ty=afft*dpix/dky;
    K=sqrt(dkx^2+dky^2);
    T=afft*dpix/K

%% DIR 2
miSIM20=SIM20;miSIM21=SIM21;miSIM22=SIM22;


    AFFT=abs(fftshift(fft2(miSIM20)))+abs(fftshift(fft2(miSIM21)))+abs(fftshift(fft2(miSIM22)));
     test=SupFun.kreis(6,size(AFFT,1),size(AFFT,2),-sind(61)*freq,-cosd(61)*freq)+SupFun.kreis(6,size(AFFT,1),size(AFFT,2),sind(61)*freq,cosd(61)*freq);
     
     
     AFFT=AFFT.*test;

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

    Tx=afft*dpix/dkx;
    Ty=afft*dpix/dky;
    K=sqrt(dkx^2+dky^2);
    T=afft*dpix/K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% start processing 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count=1; % book keeping index starting at 1

for l=1
    conW=0;  % warning flag for condition number of unmixing matrix
    


X00=SIM00;
X01=SIM01;
X02=SIM02;

X10=SIM10;
X11=SIM11;
X12=SIM12;

X20=SIM20;
X21=SIM21;
X22=SIM22;



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

[xdim,ydim]=size(X00);


%% RL deconvolution of raw data before processing
bright00=X00+X01+X02; 
if preRL
X00=deconvlucy(X00,PSF,RL1);
X01=deconvlucy(X01,PSF,RL1);
X02=deconvlucy(X02,PSF,RL1);
bright00=deconvlucy(bright00,PSF,RL1);
end
bright0=X00+X01+X02;        %% widefield image by summing images from direction 0
XF00=fftshift(fft2(X00));
XF01=fftshift(fft2(X01));
XF02=fftshift(fft2(X02));
BF00=fftshift(fft2(bright00));


%% unmixing of bands using estimated phase steps
 
 co=3
[XF00,XF01,XF02,M]=SupFun.Demod(XF00,XF01,XF02,0,phasex01,phasex01+phasex02);


disp('Demodulation of spectra finished')

%% zero pad spectra 
XF00=padarray(XF00,[round(size(XF00,2)/2),round(size(XF00,2)/2)]);
XF01=padarray(XF01,[round(size(XF01,2)/2),round(size(XF01,2)/2)]);
XF02=padarray(XF02,[round(size(XF02,2)/2),round(size(XF02,2)/2)]);
BF00=padarray(BF00,[round(size(BF00,2)/2),round(size(BF00,2)/2)]);
[xdim,ydim]=size(XF00);


%% shifting of the spectra in real space

S00f=XF00;

S01f=circshift(XF01,[-round(dky0),-round(dkx0)]);
S02f=circshift(XF02,[round(dky0),round(dkx0)]);  
  

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

mC00=SupFun.kreis(round(Kc*MiCoRa),xdim,xdim,0,0);             % Missing cone infill
mC01=SupFun.kreis(round(Kc*NotchW),xdim,xdim,-dkx0,-dky0);  % Notch filter mask 1
mC02=SupFun.kreis(round(Kc*NotchW),xdim,xdim,dkx0,dky0);    % Notch filter mask 2

Apo00=SupFun.kreis(round(Kc),xdim,xdim,0,0);       % binary apodization masks
Apo01=SupFun.kreis(round(Kc),xdim,xdim,-dkx0,-dky0); % binary apodization masks
Apo02=SupFun.kreis(round(Kc),xdim,xdim,dkx0,dky0);   % binary apodization masks

S00f=S00f.*Apo00;S01f=S01f.*Apo01;S02f=S02f.*Apo02;  %binary apodization of bands


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SIM processing for Dir1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

bright01=X00+X01+X02;

if preRL
X10=deconvlucy(X10,PSF,RL1);
X11=deconvlucy(X11,PSF,RL1);
X12=deconvlucy(X12,PSF,RL1);
bright01=deconvlucy(bright01,PSF,RL1);
end

bright1=X00+X01+X02;
XF10=fftshift(fft2(X10));
XF11=fftshift(fft2(X11));
XF12=fftshift(fft2(X12));
BF01=fftshift(fft2(bright01));

[XF10,XF12,XF11,M]=SupFun.Demod(XF10,XF11,XF12,0,phasex11,phasex12+phasex11);

disp('Demodulation of spectra finished')

%% zero pad spectra
XF10=padarray(XF10,[round(size(XF10,2)/2),round(size(XF10,2)/2)]);
XF11=padarray(XF11,[round(size(XF11,2)/2),round(size(XF11,2)/2)]);
XF12=padarray(XF12,[round(size(XF12,2)/2),round(size(XF12,2)/2)]);
BF01=padarray(BF01,[round(size(BF01,2)/2),round(size(BF01,2)/2)]);
[xdim,ydim]=size(XF00);

%% shifting of the spectra in real space

[px,py] = meshgrid(0:afft-1);

S10f=XF10;


S11f=circshift(XF11,[round(dky1),round(dkx1)]);
S12f=circshift(XF12,[-round(dky1),-round(dkx1)]);


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

mC10=SupFun.kreis(round(Kc*MiCoRa),xdim,xdim,0,0);
mC11=SupFun.kreis(round(Kc*NotchW),xdim,xdim,-dkx1,-dky1);
mC12=SupFun.kreis(round(Kc*NotchW),xdim,xdim,dkx1,dky1);

Apo10=SupFun.kreis(round(Kc),xdim,xdim,0,0);
Apo11=SupFun.kreis(round(Kc),xdim,xdim,-dkx1,-dky1);
Apo12=SupFun.kreis(round(Kc),xdim,xdim,dkx1,dky1);


S10f=S10f.*Apo10;S11f=S11f.*Apo12;S12f=S12f.*Apo11;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SIM processing of Dir2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 
bright02=X20+X21+X22;
if preRL
X20=deconvlucy(X20,PSF,RL1);
X21=deconvlucy(X21,PSF,RL1);
X22=deconvlucy(X22,PSF,RL1);
bright02=deconvlucy(bright02,PSF,RL1);
end

bright2=X20+X21+X22;
XF20=fftshift(fft2(X20));
XF21=fftshift(fft2(X21));
XF22=fftshift(fft2(X22));
BF02=fftshift(fft2(bright02));

%% unmix bands

[XF20,XF21,XF22,M]=SupFun.Demod(XF20,XF21,XF22,0,phasex21,phasex22+phasex21);


disp('Demodulation of spectra finished')

%% zero pad spectra
XF20=padarray(XF20,[round(size(XF20,2)/2),round(size(XF20,2)/2)]);
XF21=padarray(XF21,[round(size(XF21,2)/2),round(size(XF21,2)/2)]);
XF22=padarray(XF22,[round(size(XF22,2)/2),round(size(XF22,2)/2)]);
BF02=padarray(BF02,[round(size(BF02,2)/2),round(size(BF02,2)/2)]);
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

mC20=SupFun.kreis(round(Kc*MiCoRa),xdim,xdim,0,0);
mC21=SupFun.kreis(round(Kc*NotchW),xdim,xdim,-dkx2,-dky2);
mC22=SupFun.kreis(round(Kc*NotchW),xdim,xdim,dkx2,dky2);

Apo20=SupFun.kreis(round(Kc),xdim,xdim,0,0);
Apo21=SupFun.kreis(round(Kc),xdim,xdim,-dkx2,-dky2);
Apo22=SupFun.kreis(round(Kc),xdim,xdim,dkx2,dky2);


S20f=S20f.*Apo20;S21f=S21f.*Apo21;S22f=S22f.*Apo22;


%% Assemble spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



m_out=4;
mod=5;

ApoTot=Apo00+mod*Apo01+mod*Apo02+Apo10+mod*Apo11+mod*Apo12+Apo20+mod*Apo21+mod*Apo22;
indS=find(ApoTot==mod);

ApoTot(indS)=mod*m_out;

indS=find(ApoTot==2*mod);
ApoTot(indS)=ApoTot(indS)*2;
indS=find(ApoTot==2*mod+3);
ApoTot(indS)=ApoTot(indS)*1.8;

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
Recon(ind00)=(S01f(ind00)+S11f(ind00)+S21f(ind00)+S22f(ind00)+S12f(ind00)+S02f(ind00)+(S00f(ind00)+S20f(ind00)+S10f(ind00))*WFweight)*CenterWeight;
Recon(ind00)=(S01f(ind00)+S11f(ind00)+S21f(ind00)+S22f(ind00)+S12f(ind00)+S02f(ind00)+(BF00(ind00)+BF01(ind00)+BF02(ind00))*WFweight)*CenterWeight;
end

% notch filter / repace pattern peaks with zero order band information
Recon(ind01)=S00f(ind01)*Nweight;
Recon(ind02)=S00f(ind02)*Nweight;
Recon(ind11)=S00f(ind11)*Nweight;
Recon(ind12)=S00f(ind12)*Nweight;
Recon(ind21)=S00f(ind21)*Nweight;
Recon(ind22)=S00f(ind22)*Nweight;
%%




Recon2=real(ifft2(ifftshift(Recon)));  

Recon2=deconvlucy(Recon2,PSFsim,RL1);

Recon3=fftshift(fft2(Recon2));

T=SupFun.Triangle(scaleT*round(K+Kc),xdim,xdim,0,0);
if apod
Recon3=Recon3.*T;
end
Recon3=real(ifft2(ifftshift(Recon3))); 
 
if plo
RLRecon=(fftshift(fft2(Recon3)));
figure; imshow(log(abs(Recon)),[])
figure; imshow(log(abs(RLRecon)),[])
figure;imshow(Recon3,[])
figure;imshow(imresize(WFSum,2),[])
end

bl=round(size(Recon3,1)/4);
if plo
figure; subplot(1,2,1);imshow(Recon3(bl:end-bl,bl:end-bl),[]);

tempW=imresize(WFSum,2);
subplot(1,2,2);imshow(tempW(bl:end-bl,bl:end-bl),[])

end

%% get the brightfield spectra

bright=XF00+XF01+XF02+XF10+XF11+XF12+XF20+XF21+XF22;


%% inverse FFT of reconstructed data/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/OPMSIM/projection/Reto/U2OS/MIC60/220504/Cell5_BJtest/
bright=real(ifft2(ifftshift(bright)));


disp(['plane ' num2str(count) ' processed'])
count=count+1

end

%% saving of data
if saving
   
RecStack=Recon3;

RecStack=RecStack/max(RecStack(:));

WFSum=WFSum/max(WFSum(:));

RecStack=RecStack*20000;

WFSum=WFSum*20000;


imwrite(uint16(WFSum),[recFolder 'WidefieldstackSum-' num2str(n) '.tif']);

imwrite(uint16(RecStack),[recFolderM '/SIM-' num2str(n) '.tif']);


end


end

