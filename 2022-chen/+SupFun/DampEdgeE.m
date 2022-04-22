%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Program to damp Gibb's Oscillation in FFTs
%   mixing one side into the other, with a gaussian weight
%
%   by Reto Fiolka, 7.11.2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

function c=DampEdgeC(a,w)%Bild einlesen
a=double(a);
b=double(a);
%w=40;  %weite des zu mixenden Randes
xdim=size(a,2);
%Konstruktion der Gewichtungsfunktion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=1:2*w;
we=exp(-(x-w).^2/(round(w/2)^2))/2+0.5;
weL=we(1:w);
weR=we(w+1:end);
weLI=1-weL; %inverse Gewichtung
weRI=1-weR;


%a=a(1:1024,1:1024);

bordL=a(:,end-w:end); %dieser Rand wird auf der anderen Seite dann reingemixt
bordR=a(:,1:w);
bordR=a(:,1:w);
bordL=a(:,end-w+1:end);
bordU=a(1:w,:);
bordO=a(end-w+1:end,:);
bordU=flipud(bordU); %spiegeln
bordO=flipud(bordO);
bordL=fliplr(bordL);
bordR=fliplr(bordR);

%Gewichtung der R�nder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bordL=bordL.*repmat(weLI,[xdim,1]);
a(:,end-w+1:end)=a(:,end-w+1:end).*repmat(weR,[xdim,1]);

bordR=bordR.*repmat(weRI,[xdim,1]);
a(:,1:w)=a(:,1:w).*repmat(weL,[xdim,1]);

bordO=bordO.*repmat(weLI',[1,xdim]);
a(end-w+1:end,:)=a(end-w+1:end,:).*repmat(weR',[1,xdim]);

bordU=bordU.*repmat(weRI',[1,xdim]);
a(1:w,:)=a(1:w,:).*repmat(weL',[1,xdim]);

%reinmixen der R�nder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 a(:,end-w+1:end)=(a(:,end-w+1:end)+bordR);
 a(:,1:w)=(a(:,1:w)+bordL);

 a(end-w+1:end,:)=(a(end-w+1:end,:)+bordU);
 a(1:w,:)=(a(1:w,:)+bordO);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Im �berlapp nehme ich die Original Bildinformation
a(1:w,1:w)=b(1:w,1:w);
a(end-w:end,1:w)=b(end-w:end,1:w);
a(1:w,end-w:end)=b(1:w,end-w:end);
a(end-w:end,end-w:end)=b(end-w:end,end-w:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=a;