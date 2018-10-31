function [M] = unknownToNan(M,numKnownCols)
[numRow,numCol] = size(M);
M(repmat([1:numCol],numRow,1) > repmat(numKnownCols,1,numCol)) = nan;