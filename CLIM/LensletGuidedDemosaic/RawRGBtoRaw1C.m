function [ Raw1C, ColorPosition] = RawRGBtoRaw1C(RawRGB, DecodeOptions)
%RAWRGBTORAW1C

[M_raw, N_raw, ~] = size(RawRGB);
Raw1C = zeros(M_raw,N_raw, DecodeOptions.Precision);
if strcmp(DecodeOptions.DemosaicOrder, 'bggr')
    Raw1C(1:2:M_raw,1:2:N_raw)=RawRGB(1:2:M_raw,1:2:N_raw,3); %position of blue pixels
    Raw1C(1:2:M_raw,2:2:N_raw)=RawRGB(1:2:M_raw,2:2:N_raw,2); %position of green pixels
    Raw1C(2:2:M_raw,2:2:N_raw)=RawRGB(2:2:M_raw,2:2:N_raw,1); %position of red pixels
    Raw1C(2:2:M_raw,1:2:N_raw)=RawRGB(2:2:M_raw,1:2:N_raw,2); %position of green pixels
    if nargout > 1
        ColorPosition = zeros(size(RawRGB), 'single');
        ColorPosition(1:2:M_raw,1:2:N_raw,3)= 1; %position of blue pixels
        ColorPosition(1:2:M_raw,2:2:N_raw,2)= 1; %position of green pixels
        ColorPosition(2:2:M_raw,2:2:N_raw,1)= 1; %position of red pixels
        ColorPosition(2:2:M_raw,1:2:N_raw,2)= 1; %position of green pixels
    end
elseif strcmp(DecodeOptions.DemosaicOrder, 'rggb')
    Raw1C(1:2:M_raw,1:2:N_raw)=RawRGB(1:2:M_raw,1:2:N_raw,1); %position of red pixels
    Raw1C(1:2:M_raw,2:2:N_raw)=RawRGB(1:2:M_raw,2:2:N_raw,2); %position of green pixels
    Raw1C(2:2:M_raw,2:2:N_raw)=RawRGB(2:2:M_raw,2:2:N_raw,3); %position of blue pixels
    Raw1C(2:2:M_raw,1:2:N_raw)=RawRGB(2:2:M_raw,1:2:N_raw,2); %position of green pixels
    if nargout > 1
        ColorPosition = zeros(size(RawRGB), 'single');
        ColorPosition(1:2:M_raw,1:2:N_raw,1)= 1; %position of red pixels
        ColorPosition(1:2:M_raw,2:2:N_raw,2)= 1; %position of green pixels
        ColorPosition(2:2:M_raw,2:2:N_raw,3)= 1; %position of blue pixels
        ColorPosition(2:2:M_raw,1:2:N_raw,2)= 1; %position of green pixels
    end
elseif strcmp(DecodeOptions.DemosaicOrder, 'gbrg')
    Raw1C(1:2:M_raw,1:2:N_raw)=RawRGB(1:2:M_raw,1:2:N_raw,2); %position of green pixels
    Raw1C(1:2:M_raw,2:2:N_raw)=RawRGB(1:2:M_raw,2:2:N_raw,3); %position of blue pixels
    Raw1C(2:2:M_raw,2:2:N_raw)=RawRGB(2:2:M_raw,2:2:N_raw,2); %position of green pixels
    Raw1C(2:2:M_raw,1:2:N_raw)=RawRGB(2:2:M_raw,1:2:N_raw,1); %position of red pixels
    if nargout > 1
        ColorPosition = zeros(size(RawRGB), 'single');
        ColorPosition(1:2:M_raw,1:2:N_raw,2)= 1; %position of green pixels
        ColorPosition(1:2:M_raw,2:2:N_raw,3)= 1; %position of blue pixels
        ColorPosition(2:2:M_raw,2:2:N_raw,2)= 1; %position of green pixels
        ColorPosition(2:2:M_raw,1:2:N_raw,1)= 1; %position of red pixels
    end
elseif strcmp(DecodeOptions.DemosaicOrder, 'grbg')
    Raw1C(1:2:M_raw,1:2:N_raw)=RawRGB(1:2:M_raw,1:2:N_raw,2); %position of green pixels
    Raw1C(1:2:M_raw,2:2:N_raw)=RawRGB(1:2:M_raw,2:2:N_raw,1); %position of red pixels
    Raw1C(2:2:M_raw,2:2:N_raw)=RawRGB(2:2:M_raw,2:2:N_raw,2); %position of green pixels
    Raw1C(2:2:M_raw,1:2:N_raw)=RawRGB(2:2:M_raw,1:2:N_raw,3); %position of blue pixels
    if nargout > 1
        ColorPosition = zeros(size(RawRGB), 'single');
        ColorPosition(1:2:M_raw,1:2:N_raw,2)= 1; %position of green pixels
        ColorPosition(1:2:M_raw,2:2:N_raw,1)= 1; %position of red pixels
        ColorPosition(2:2:M_raw,2:2:N_raw,2)= 1; %position of green pixels
        ColorPosition(2:2:M_raw,1:2:N_raw,2)= 1; %position of blue pixels
    end
else
    fprintf('Error : BayerPattern not valid\n');
end

end

