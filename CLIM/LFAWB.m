function [out] = LFAWB(im_orig,LF,option,T,maxIter)
%robustAWB(filename,out,option,catType,T,maxIter,plot)
% Performs robust auto white-balancing by estimating gray pixels based on
% their deviation in YUV space then applying an iterative correction via
% CAT or directly adjusting the R and B channels.
% Set option = 'RB gain' or 'cat' for the correction method.
% Set T higher for a larger threshold of deviation to consider off-gray.

% tutorialinit

%% import image
% im_orig = cbimread(filename);
% figure
% imshow(im_orig)
% title('Original Image')

%% various constants
xyz_D65 = [95.04; 100; 108.88]; %http://en.wikipedia.org/wiki/D65, normalized Y = 100

sRGBtoXYZ =   [0.4124564  0.3575761  0.1804375; ...
               0.2126729  0.7151522  0.0721750; ...
               0.0193339  0.1191920  0.9503041];

%% robust auto white balance - YUV grays
% option = 'RB gain'; % adjustments are done as a gain to only R or B channelas descripted in Huo
% option = 'cat';
% catType = 'CAT02';

% maxIter = 1000;

% T = 0.3; %this should be tweakable slider
% T = 0.1;
% gain adjustment parameters, can probably be optimized via more control system analysis
u = .01; %gain step size
a = .8; %double step threshold
b = .001; %convergence threshold

% rgb to yuv
xfm =   [0.299 0.587 0.114; ...
        -0.299 -0.587 0.886; ...
        0.701 -0.587 -0.114];
    
inv_xfm = inv(xfm);

LF = LF';
im_orig = im_orig';

% imRGB_orig = cbreshape(im_orig)*255;
imRGB_orig = im_orig;
imRGB = imRGB_orig;
gain = [1 1 1];
catm_acc = eye(3);

U_avg = [];
V_avg = [];
totGray = [];
for iter = 1:maxIter
    im = xfm*imRGB; %convert to YUV
    % find gray chromaticity
    % (|U|+|V|)/Y
    F = ( abs(im(2,:)) + abs(im(3,:)) )./im(1,:);
    
%     figure
%     imshow(cbunshape(imRGB/255.*repmat(F<T,3,1),size(im_orig)))
%     % imshow(cbunshape(inv_xfm)*im+[F>T;zeros(2,length(F))],size(im_orig))) %rubylith overlays
%     title('Gray Pixels Found by Robust AWB'),xlabel('Pixels in their original color are considered gray')
    
    totGray = [totGray sum(F<T)];
    if totGray(end) == 0
        disp('No valid gray pixels found.')
        break
    end
    grays = im(:,F<T);
    U_bar = mean(grays(2,:));
    V_bar = mean(grays(3,:));
    
    U_avg = [U_avg U_bar];
    V_avg = [V_avg V_bar];
    
    if strcmpi(option,'cat')
        if max(abs([U_bar V_bar])) < b
            disp(['Converged. U_bar and V_bar < ' num2str(b) ' in magnitude.'])
            break
        elseif iter >= 2 && norm([U_avg(end)-U_avg(end-1) V_avg(end)-V_avg(end-1)]) < 10^-6
            disp(['U_bar and V_bar are no longer improving.'])
            break
        end
        rgbEst = inv_xfm*[100;U_bar;V_bar]; %convert the average gray from YUV to RGB
        xyEst = LFXYZ2xy(sRGBtoXYZ*rgbEst); %calculate xy chromaticity
        xyzEst = LFxy2XYZ(xyEst,100); %normalize Y to 100 so D65 luminance comparable
        catm = LFcbCAT(xyzEst,xyz_D65);
        imRGB = catm*imRGB;
%         LF = catm * LF;
        catm_acc = catm * catm_acc;
    else
        if abs(U_bar) > abs(V_bar) % U > V; blue needs adjustment
            err = U_bar;
            ch = 3; %blue channel
        else
            err = V_bar;
            ch = 1; %red channel
        end
        if abs(err) >= a
            delta = 2*sign(err)*u; %accelerate gain adjustment if far off
        elseif abs(err) < b %converged
            delta = 0;
            disp(['Converged. U_bar and V_bar < ' num2str(b) ' in magnitude.'])
            break
        else
            delta = err*u;
        end
%         [err,delta]
        gain(ch) = gain(ch)-delta; %negative fdbk loop
        imRGB = diag(gain)*imRGB_orig;
    end
end

out = diag(gain)*LF;
out = (catm_acc * out)';

end

function outMat = LFcbCAT(xyz_est,xyz_target)
%outMat = cbCAT(xyz_est,xyz_target,type)
% Chromatic adaptation transform via von Kries's method.
%  type     chooses the LMS-like space to apply scaling in, valid options:
%           'vonKries', 'bradford', 'sharp', 'cmccat2000', 'cat02', 'xyz'
% See http://www.brucelindbloom.com/index.html?Eqn_ChromAdapt.html

xyz_est = LFmakecol(xyz_est);
xyz_target = LFmakecol(xyz_target);

% the following are mostly taken from S. Bianco. "Two New von Kries Based
% Chromatic Adapatation Transforms Found by Numerical Optimization."

%Hunt-Pointer-Estevez normalized to D65
xfm = [.40024 .7076 -.08081; -.2263 1.16532 0.0457; 0 0 .91822];

    
gain = (xfm*xyz_target)./(xfm*xyz_est);
outMat = xfm\(diag(gain)*xfm);
sRGBtoXYZ =   [0.4124564  0.3575761  0.1804375; ...
               0.2126729  0.7151522  0.0721750; ...
               0.0193339  0.1191920  0.9503041];
outMat = inv(sRGBtoXYZ)*outMat*sRGBtoXYZ;
end

function x = LFmakecol(x)
%x = makecol(x)
% returns x as a column vector

s = size(x);
if (length(s) == 2) && (s(1) < s(2))
    x = x.';
end

end

function xy = LFXYZ2xy(xyz)
%xy = XYZ2xy(xyz)
% Converts CIE XYZ to xy chromaticity.

X = xyz(1);
Y = xyz(2);
s = sum(xyz);
xy = [X/s; Y/s];

end

function xyz = LFxy2XYZ(xy,Y)
%xyz = xy2XYZ(xy,Y)
% Converts xyY chromaticity to CIE XYZ.
x = xy(1); y = xy(2);
xyz = [Y/y*x; Y; Y/y*(1-x-y)];
end
