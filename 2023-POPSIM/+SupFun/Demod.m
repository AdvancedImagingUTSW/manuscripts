function [S1,S2,S3,M]=Demod(X0,X1,X2,phasex0,phasex1,phasex2)

%This function solves the linear equation system for HELM
%
% Reto Fiolka, 17.07.2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=[X0; X1; X2]; 
i=sqrt(-1);
m=2; 
M=[m exp(i*phasex0) exp(-i*phasex0);
   m exp(i*phasex1) exp(-i*phasex1); 
   m exp(i*phasex2) exp(-i*phasex2)]; Kondition=cond(M)
Minv=inv(M);

S=zeros(size(X));
[a1,b1]=size(X0);
for i=1:a1
    temp=Minv*([X(i,1:b1); X(i+a1,1:b1); X(i+2*a1,1:b1)]);
    S(i,:)=temp(1,:);S(a1+i,:)=temp(2,:); S(2*a1+i,:)=temp(3,:); 
end
clear('X','X0','X1','X2')
S1=S(1:a1,1:b1);        %original spectra
S2=S(1+a1:2*a1,1:b1);   %shifted spectra (+u)
S3=S(1+2*a1:3*a1,1:b1); %shifted spectra (-u)

end
