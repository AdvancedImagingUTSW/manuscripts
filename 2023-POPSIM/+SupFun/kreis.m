
function cut=kreis(r,xdim,ydim,midy,midx)
cut=zeros(xdim,ydim);
for i = 1:xdim
    for j = 1:ydim

        x=i-xdim/2-midx;  y=j-ydim/2-midy;

        %disk
        if((x^2+y^2)<r^2)
            cut(i,j)=cut(i,j)+1;
        end

    end
end

end


