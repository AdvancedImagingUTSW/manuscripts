
function T=Triangle3D(r,xdim,ydim,zdim,midy,midx,midz)
T=zeros(xdim,ydim,zdim);
for i = 1:xdim
    for j = 1:ydim
        for k = 1:zdim

        x=i-xdim/2-midx;  
        y=j-ydim/2-midy;
        z=k-round(zdim/2)-midz;
          R=sqrt((x^2+y^2+z^2));
        %disk
        if((x^2+y^2+(z*xdim/zdim*2.5)^2)<r^2)
            T(i,j,k)=1-R/r;
            
        end
    end
    end
 
end


