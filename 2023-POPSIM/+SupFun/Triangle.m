
function T=Triangle(r,xdim,ydim,midy,midx)
T=zeros(xdim,ydim);
for i = 1:xdim
    for j = 1:ydim

        x=i-xdim/2-midx;  y=j-ydim/2-midy;
          R=sqrt((x^2+y^2));
        %disk
        if((x^2+y^2)<r^2)
            T(i,j)=1-R/r;
        end

    end
end

end


