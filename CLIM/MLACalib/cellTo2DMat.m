% Utility function used for mla_calib metadata reading.
% Convert nx1 cells of row vectors into a 2D matrix.
%Inputs : 
% - Cell : input nx1 cell of row vectors.
% - matsize : column and row length of the matrix given as a vector.
function Mat = cellTo2DMat(Cell,matsize)
Mat = zeros(matsize);

for i=1:matsize(1)
    minSize = min(matsize(2),length(Cell{i}));
    Mat(i,1:minSize) = Cell{i}(1:minSize);
end