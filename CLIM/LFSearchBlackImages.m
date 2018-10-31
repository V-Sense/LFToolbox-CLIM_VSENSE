% LFSearchBlackImages - Searches Black images corresponding to a White image file name(including full path and extension)
%
% Usage: 
%     [ImagesDir, BlackRawNames] = LFSearchBlackImages(WhiteRawFname,BlackImageIndices);
% 
% Inputs :
% 
%    WhiteRawFname : Full name of a raw white image.
%
%    BlackImageIndices : vector containing the indices of the black images
%    in the White raw Image Folder. (Black images are assumed to be
%    in the same folder as the White images).
% 
% Outputs : 
% 
%     ImagesDir : White (and Black) image folder.
%     
%     BlackRawNames : cell containing the filenames of all the black
%     images found in the folder corresponding to the indices given in
%     BlackImageIndices.

function [ImagesDir, BlackRawNames] = LFSearchBlackImages(WhiteRawFname, BlackImageIndices)

[ImagesDir, WhiteImageName, ImageExt] = fileparts(WhiteRawFname);

%Search for format of the white image file name
[DigitsStart, DigitsEnd] = regexp(WhiteImageName,'_\d+','start','end');
DigitsStart = DigitsStart(end);
DigitsEnd = DigitsEnd(end);

%Retrieve black images corresponding to the same format and black image indices
fdir = dir(ImagesDir);
NumExpr = sprintf('%i|',floor(BlackImageIndices));
NumExpr = ['(' NumExpr(1:end-1) ')'];

BlackRawNames=regexp({fdir.name},['^' WhiteImageName(1:DigitsStart) '0*' NumExpr WhiteImageName(DigitsEnd+1:end) ImageExt '$'],'match');
BlackRawNames=[BlackRawNames{:}];