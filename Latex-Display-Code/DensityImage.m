folderPath = 'C:\Users\rusco\Downloads\Density';
imageFiles = dir(fullfile(folderPath, '*.png'));
nImages = length(imageFiles);

firstImage = imread(fullfile(folderPath, imageFiles(1).name));
if size(firstImage, 3) > 1
    firstImage = rgb2gray(firstImage);
end
imageSize = size(firstImage);
sumImage = zeros(imageSize);

% loop images
for k = 1:nImages
    currentImage = imread(fullfile(folderPath, imageFiles(k).name));
    
    % convert to grayscale if necessary
    if size(currentImage, 3) > 1
        currentImage = rgb2gray(currentImage);
    end
    
    binaryImage = imbinarize(currentImage);
    
    % sum the images together
    sumImage = sumImage + double(binaryImage);
end

% compute density, ie = 0 if no images cross the point, = 1 if all do
densityImage = sumImage / nImages;

% plot the density map
figure;
imagesc(densityImage);
colormap('jet'); % Use a colormap of your choice
colorbar;
axis image off;
title('Graded Lattice Density');

densityImage((0<densityImage) & (densityImage<0.5)) = 0.5;

% plot adjusted density map
figure;
imagesc(densityImage);
colormap('jet'); % Use a colormap of your choice
colorbar;
axis image off;
title('Adjusted Graded Lattice Density');