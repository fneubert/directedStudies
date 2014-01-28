function [U, S, D] = test()
%% Test client for Stochastic PCA on FCVG Data

% Load data mean size
load('imgDataMeanSize.mat');

% Load n sample names in to string array
samples = readStringStruct();

% Set dimensionality parameters
c = 15;
meanSize = round(meanSize./c);
d = prod([meanSize 3]);
k = 30;

% Perform sPCA
alignmentIO = 1; % With/without alignment
[U, S, D] = incrementalPCA2(samples, d, k, meanSize, alignmentIO);

% Visualize the top K eigenvectors found
displayData(U', meanSize);

end
