function meanSize = getMeanSize()

filePath = '/Users/Fariz/Documents/Kurser/ECP/directedStudies/fgvcData/data/images/';
files    = dir(fullfile(filePath, '*.jpg'));

n = numel(files);
dataSizes = zeros(n, 2);

for m = 1:n
    fileDir = strcat(filePath, files(m).name);
    img = imfinfo(fileDir, 'jpg');
    dataSizes(m, :) = [img.Width img.Height];
    fprintf('Iteration %5d: \n', m);
end

meanSize = mean(dataSizes);

end
