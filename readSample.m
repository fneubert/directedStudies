function x = readSample(sample, meanSize)

% Image filepath
filePath = '/Users/Fariz/Documents/Kurser/ECP/directedStudies/fgvcData/data/images/';

% Read, resize and vectorize image sample
fileDir = strcat(filePath, sample);
img = imread(fileDir, 'jpg');

% Resize images to mean size
imgResize = imresize(img, meanSize, 'nearest');

% Turn into (double type) row vector
x = double(imgResize);
x = x(:);
%x = normalizeData(x);
