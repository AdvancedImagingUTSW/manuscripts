function tiffWrite(imData,PSFname)

[nx, ny, nz]= size(imData);
imgType= class(imData);
tagstruct.Photometric= Tiff.Photometric.MinIsBlack;
tagstruct.ImageLength = nx;
tagstruct.ImageWidth = ny;
tagstruct.PlanarConfiguration= Tiff.PlanarConfiguration.Chunky;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.BitsPerSample= 16;
tiffFile=Tiff(PSFname, 'w');

for iz=1:nz
    tiffFile.setTag(tagstruct);
    tiffFile.write(imData(:,:,iz));
    tiffFile.writeDirectory();
    
end
tiffFile.close();