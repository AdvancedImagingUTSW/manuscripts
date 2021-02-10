function [tiffImage]=tiffRead(imagePath)

% Determine Image Properties
InfoImage=imfinfo(imagePath);

% Pre-allocate Memory
tiffImage = zeros(InfoImage(1).Height,InfoImage(1).Width,length(InfoImage),'uint16');

% Iteratively Load the Image
TifLink = Tiff(imagePath, 'r');
for i=1:length(InfoImage)
    TifLink.setDirectory(i);
    tiffImage(:,:,i)=TifLink.read();
end
TifLink.close();
