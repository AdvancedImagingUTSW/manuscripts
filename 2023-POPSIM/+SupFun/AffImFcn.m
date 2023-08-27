function Cor=AffImFcn(x,Im1,Im2)

tform = affine2d([x(1) 0 0; 0 x(2) 0; 0 0 1]);
temp = imwarp(Im1,tform);
temp2=imrotate(temp,x(3));
 c = normxcorr2(Im2,temp2);
 
 Cor=1-max(c(:));
end


