% LFColourCorrect - applies a colour correction matrix, balance vector, and gamma, called by LFUtilDecodeLytroFolder
%
% Usage: 
%     LF = LFColourCorrect( LF, ColMatrix, ColBalance, Gamma )
% 
% This implementation deals with saturated input pixels by aggressively saturating output pixels.
%
% Inputs :
% 
%     LF : a light field or image to colour correct. It should be a floating point array, and
%          may be of any dimensinality so long as the last dimension has length 3. For example, a 2D
%          image of size [Nl,Nk,3], a 4D image of size [Nj,Ni,Nl,Nk,3], and a 1D list of size [N,3]
%          are all valid.
%
%    ColMatrix : a 3x3 colour conversion matrix. This can be built from the metadata provided
%                with Lytro imagery using the command:
%                     ColMatrix = reshape(cell2mat(LFMetadata.image.color.ccmRgbToSrgbArray), 3,3);
%                as demonstrated in LFUtilDecodeLytroFolder.
% 
%    ColBalance : 3-element vector containing a multiplicative colour balance.
% 
%    Gamma : rudimentary gamma correction is applied of the form LF = LF.^Gamma.
%
% Outputs : 
% 
%     LF, of the same dimensionality as the input.
% 
% 
% See also: LFHistEqualize, LFUtilDecodeLytroFolder

% Part of LF Toolbox v0.4 released 12-Feb-2015
% Copyright (c) 2013-2015 Donald G. Dansereau

% modified by Rodrigo Daudt : 22 Aug. 2016 :
%   -added Automatic white balance.

function LF = LFColourCorrect(LF, ColMatrix, ColBalance, ExposureBias, Gamma, DoAWB, noClip)

LFSize = size(LF);

% Flatten input to a flat list of RGB triplets
NDims = numel(LFSize);
c1 = ceil(LFSize(1)/2);
c2 = ceil(LFSize(2)/2);
c_slice = squeeze(LF(c1,c2,:,:,:));
c_slice = reshape(c_slice, [prod(LFSize(3:NDims-1)), 3]);
LF = reshape(LF, [prod(LFSize(1:NDims-1)), 3]);


LF = bsxfun(@times, LF, ColBalance);
LF = LF * ColMatrix;

% noClip keeps all the pixel values above the saturation level in the pipeline.
SaturationLevel = 2^(-ExposureBias);
if(~noClip)
    LF = min(SaturationLevel,max(0,LF)) ./ SaturationLevel;
else
    LF = max(0,LF) ./ SaturationLevel;
end

% Apply gamma
if(isa(Gamma,'char') && strcmp('sRGB',Gamma))
    LF = sRGB_Gamma(LF);
elseif(isnumeric(Gamma))
    LF = LF .^ Gamma;
else
    warning('unrecognized format for gamma parameter -> default sRGB correction applied.');
    LF = sRGB_Gamma(LF);
end

% Apply AWB
if(DoAWB)
	LF = LFAWB(c_slice,LF,'cat',0.35,10000);
	%LF = LFAWB(c_slice,LF,'RB gain',0.35,10000);
end

% Unflatten result
LF = reshape(LF, [LFSize(1:NDims-1),3]);

end



function Iout = sRGB_Gamma(Iin)
    Mask = Iin < 0.0031308;
    Iout = max(0,Iin*12.92) .* Mask + (1.055*realpow(max(Iin,0),1/2.4) - 0.055) .* (~Mask);
end