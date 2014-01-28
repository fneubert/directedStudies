function xt = alignImage(xt, U, S, meanSize, K)
% Determines whether the flipped or non-flipped image received the
% highest scores and returns the one with the highest
[xt1,score1] = translateImage(xt, U, S, meanSize, K);
[xt,score2] = translateImage(flipImage(xt,meanSize), U, S, meanSize, K);

if(score1>score2)
    xt = xt1;
end
end

function xt = flipImage(xt, meanSize)
% Flips images from left to right
xt = reshape(xt, [meanSize 3]);
for i=1:3
    xt(:,:,i) = fliplr(xt(:,:,i));
end
end

    
function [xt,score] = translateImage(xt, U, S, meanSize, K)
% Calculates the projections of all allowed translations of image
% xt onto the current eigenspace spanned by the columns of U

m = meanSize(1); n = meanSize(2);
X = reshape(xt, [meanSize 3]);
C = 0; %zeros(2*m - 1, 2*n - 1);

ypad = round(m/2);
xpad = round(n/2);

for i = 1:K
    Ui = reshape(U(:, i), [meanSize 3]);
    for j = 1:3
        Xj = X(:, :, j);
        Uij = Ui(:, :, j);
        %C = C + (xcorr2(Xj, Uij).^2).*S(i, i);
        C = C + (xcorr2(Xj, padarray(Uij,[ypad,xpad],'circular')).^2);
    end
end

C = C(ypad+1:end-ypad,xpad+1:end-xpad);

score = max(C(:));

[k, l] = find(C == max(C(:)));

k = k(1) - m;
l = l(1) - n;

xt = zeros(m,n,3);
xt(max(1, k+1):min(m, m+k),max(1, l+1):min(n, n+l),:) = ...
    X(max(1, 1-k):min(m, m-k),max(1, 1-l):min(n, n-l),:);

xt = xt(:);

end
