function [ ColorMatrix] = GetBayerMatrix(M,N, BayerPattern )
%GETBAYERMATRIX Summary of this function goes here
%   Detailed explanation goes here

ColorMatrix = zeros(M,N, 'single');
if strcmp(BayerPattern, 'bggr')        
        ColorMatrix(1:2:M,1:2:N)= 3; %position of blue pixels
        ColorMatrix(1:2:M,2:2:N)= 2; %position of green pixels
        ColorMatrix(2:2:M,2:2:N)= 1; %position of red pixels
        ColorMatrix(2:2:M,1:2:N)= 2; %position of green pixels

elseif strcmp(BayerPattern, 'rggb')
        ColorMatrix(1:2:M,1:2:N)= 1; %position of red pixels
        ColorMatrix(1:2:M,2:2:N)= 2; %position of green pixels
        ColorMatrix(2:2:M,2:2:N)= 3; %position of blue pixels
        ColorMatrix(2:2:M,1:2:N)= 2; %position of green pixels

elseif strcmp(BayerPattern, 'gbrg')
        ColorMatrix(1:2:M,1:2:N)= 2; %position of green pixels
        ColorMatrix(1:2:M,2:2:N)= 3; %position of blue pixels
        ColorMatrix(2:2:M,2:2:N)= 2; %position of green pixels
        ColorMatrix(2:2:M,1:2:N)= 1; %position of red pixels

elseif strcmp(BayerPattern, 'grbg')
        ColorMatrix(1:2:M,1:2:N)= 2; %position of green pixels
        ColorMatrix(1:2:M,2:2:N)= 1; %position of red pixels
        ColorMatrix(2:2:M,2:2:N)= 2; %position of green pixels
        ColorMatrix(2:2:M,1:2:N)= 3; %position of blue pixels
else
    fprintf('Error : BayerPattern not valid\n');
end

end

