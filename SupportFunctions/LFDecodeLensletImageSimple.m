% LFDecodeLensletImageSimple - decodes a 2D lenslet image into a 4D light field, called by LFUtilDecodeLytroFolder
%
% Usage:
%
%   [LF, LFWeight, DecodeOptions, DebayerLensletImage, CorrectedLensletImage] = ...
%      LFDecodeLensletImageSimple( LensletImage, WhiteImage, LensletGridModel, DecodeOptions )
%
% This function follows the simple decode process described in:
% D. G. Dansereau, O. Pizarro, and S. B. Williams, "Decoding, calibration and rectification for
% lenslet-based plenoptic cameras," in Computer Vision and Pattern Recognition (CVPR), IEEE
% Conference on. IEEE, Jun 2013.
%
% This involves demosaicing, devignetting, transforming and slicing the input lenslet image to
% yield a 4D structure. More sophisticated approaches exist which combine steps into joint
% solutions, and they generally yield superior results, particularly near the edges of lenslets.
% The approach taken here was chosen for its simplicity and flexibility.
%
% Input LensletImage is a raw lenslet iamge as found within Lytro's LFP picture files. Though
% this function is written with Lytro imagery in mind, it should be possible to adapt it for use
% with other lenslet-based cameras. LensletImage is ideally of type uint16.
%
% Input WhiteImage is a white image as found within Lytro's calibration data. A white image is an
% image taken through a diffuser, and is useful for removing vignetting (darkening near edges of
% images) and for locating lenslet centers. The white image should be taken under the same zoom and
% focus settings as the light field. In the case of Lytro imagery, the helper functions
% LFUtilProcessWhiteImages and LFSelectFromDatabase are provided to assist in forming a list of
% available white images, and selecting the appropriate image based on zoom and focus settings.
%
% Input LensletGridModel is a model of the lenslet grid, as estimated, for example, using
% LFBuildLensletGridModel -- see that function for more on the required structure.
%
% Optional Input DecodeOptions is a structure containing:
%          [Optional] ResampMethod : 'fast'(default) or 'triangulation', the latter is slower
%                                    Or 'barycentric' : result in higher resolution (multiplied 
%                                    by 3*sqrt(3)/2 compared to the 2 other methods)
%          [Optional]   Precision : 'single'(default) or 'double'
%          [Optional] LevelLimits : a two-element vector defining the black and white levels
%
% Output LF is a 5D array of size [Nj,Ni,Nl,Nk,3]. See [1] and the documentation accompanying this
% toolbox for a brief description of the light field structure.
%
% Optional output LFWeight is of size [Nj,Ni,Nl,Nk]. LFWeight contains a confidence measure
% suitable for use in filtering applications that accept a weight parameter. This parameter is kept
% separate from LF rather than building in as a fourth channel in order to allow an optimization:
% when the parameter is not requested, it is not computed, saving significant processing time and
% memory.
%
% Optional output DebayerLensletImage is useful for inspecting intermediary results. The debayered
% lenslet image is the result of devignetting and debayering, with no further processing. Omitting
% this output variable saves memory.
%
% Optional output CorrectedLensletImage is useful for inspecting intermediary results. The
% corrected lenslet image has been rotated and scaled such that lenslet centers lie on straight lines, and
% every lenslet center lies on an integer pixel spacing. See [1] for more detail. Omitting
% this output variable saves memory.
%
% See also:  LFLytroDecodeImage, LFUtilDecodeLytroFolder

% Part of LF Toolbox v0.4 released 12-Feb-2015
% Copyright (c) 2013-2015 Donald G. Dansereau

% modified by Mikael Le Pendu
%   - Added barycentric interpolation.
%   - Added Highlights Processing + WhiteImage normalization + Early White Balance.
% modified by Pierre David
%   - Added white lenslet image guided demosaicing and interpolation

function [LF, LFWeight, DecodeOptions, DebayerLensletImage, CorrectedLensletImage] = ...
    LFDecodeLensletImageSimple(LensletImage, WhiteImage, HotPixels, LensletGridModel, DecodeOptions)

%---Defaults---
DecodeOptions = LFDefaultField( 'DecodeOptions', 'LevelLimits', [min(WhiteImage(:)), max(WhiteImage(:))] );
DecodeOptions = LFDefaultField( 'DecodeOptions', 'ResampMethod', 'fast' ); %'fast', 'triangulation', 'barycentric'
DecodeOptions = LFDefaultField( 'DecodeOptions', 'Precision', 'single' );
DecodeOptions = LFDefaultField( 'DecodeOptions', 'DoDehex', true );
DecodeOptions = LFDefaultField( 'DecodeOptions', 'DoSquareST', true );
DecodeOptions = LFDefaultField( 'DecodeOptions', 'WeightedDemosaic', false ); %P. DAVID: By default WeightedDemosaic is inactive (may increase colour noise)
DecodeOptions = LFDefaultField( 'DecodeOptions', 'WeightedInterp', true ); %By default, use microlens based weighting for interpolations in the rotation/scaling.

%---Rescale image values, remove black level---
DecodeOptions.LevelLimits = cast(DecodeOptions.LevelLimits, DecodeOptions.Precision);
BlackLevel = DecodeOptions.LevelLimits(1);
WhiteLevel = DecodeOptions.LevelLimits(2);
WhiteImage = cast(WhiteImage, DecodeOptions.Precision);
WhiteImage = (WhiteImage - BlackLevel) ./ (WhiteLevel - BlackLevel);
LensletImage = cast(LensletImage, DecodeOptions.Precision);
LensletImage = (LensletImage - BlackLevel) ./ (WhiteLevel - BlackLevel);

%Define variable for bayer pattern
if(strcmp(DecodeOptions.DemosaicOrder,'grbg'))
    G1stY=1;G1stX=1;
    RstY=1;RstX=2;
    BstY=2;BstX=1;
    G2stY=2;G2stX=2;
elseif(strcmp(DecodeOptions.DemosaicOrder,'bggr'))
    BstY=1;BstX=1;
    G1stY=1;G1stX=2;
    G2stY=2;G2stX=1;
    RstY=2;RstX=2;
else
    error('unrecognized Demosaic Order');
end


%---White Image normalization for devignetting---
  %Color response normalization:
WhiteImage(RstY:2:end,RstX:2:end) = WhiteImage(RstY:2:end,RstX:2:end) * DecodeOptions.SensorNormalizeRBGains(1);%Red component
WhiteImage(BstY:2:end,BstX:2:end) = WhiteImage(BstY:2:end,BstX:2:end) * DecodeOptions.SensorNormalizeRBGains(2);%Blue component
  %Global normalization (=> value 1 at microlens centers which are not supposed to have vignetting).
WhiteImage = WhiteImage./prctile(WhiteImage(:),99.9);

%---Devignetting---
LensletImage = LensletImage ./ WhiteImage;


%Separate R,G and B pixels for the next steps (Correction of highlights + White Balance)
if(DecodeOptions.CorrectSaturated || DecodeOptions.EarlyWhiteBalance)
    R=LensletImage(RstY:2:end,RstX:2:end);
    G1=LensletImage(G1stY:2:end,G1stX:2:end);
    G2=LensletImage(G2stY:2:end,G2stX:2:end);
    B=LensletImage(BstY:2:end,BstX:2:end);
end

%--- Highlights Processing ---
if(DecodeOptions.CorrectSaturated)
    fact=max(DecodeOptions.ColourBalance);%RGB dependent factor to force pixels saturated on all components to be white after the White Balance step.
    ThSatCol=.99;%saturation threshold
    
    %Detect pixels saturated on the sensor (i.e. close to 1 before division by White image).
    RSat = find(R>ThSatCol./WhiteImage(RstY:2:end,RstX:2:end));
    G1Sat = find(G1>ThSatCol./WhiteImage(G1stY:2:end,G1stX:2:end));
    G2Sat = find(G2>ThSatCol./WhiteImage(G2stY:2:end,G2stX:2:end));
    BSat =  find(B>ThSatCol./WhiteImage(BstY:2:end,BstX:2:end));

    %Determine the component with lowest value (the last one to saturate).
    [MinSat,MinId] = min(cat(3,R,G1,G2,B),[],3);
    %Compute the value of the lowest component after White Balance.
    ColMult = DecodeOptions.ColourBalance([1,2,2,3]);
    MinSatBal = MinSat.*ColMult(MinId);
    %Compute weight indicating the amount of saturation on the sensor (i.e. before devignetting) . If all components are saturated on the sensor, the weight is 1.
    MinSat = min(MinSat.*(WhiteImage(G1stY:2:end,G1stX:2:end)+WhiteImage(G2stY:2:end,G2stX:2:end)+WhiteImage(RstY:2:end,RstX:2:end)+WhiteImage(BstY:2:end,BstX:2:end))/4,1).^2;
    clear MinId

    %Final formula combining two behaviours : only left term when at least one component is far from saturation on the sensor/ only right term when All R,G1,G2 and B pixels are saturated on the sensor.
    %In both cases inverse white balance is applied (division by WB coeff) because the White Balance will be applied later on the whole image.
    R(RSat) =   max( ( MinSatBal(RSat) .*(1-MinSat(RSat))  +  R(RSat)*fact .* MinSat(RSat) )/DecodeOptions.ColourBalance(1) ,  R(RSat));
    G1(G1Sat) = max( ( MinSatBal(G1Sat) .*(1-MinSat(G1Sat)) + G1(G1Sat)*fact .* MinSat(G1Sat) )/DecodeOptions.ColourBalance(2) , G1(G1Sat)) ;
    G2(G2Sat) = max( ( MinSatBal(G2Sat) .*(1-MinSat(G2Sat)) + G2(G2Sat)*fact .* MinSat(G2Sat) )/DecodeOptions.ColourBalance(2) , G2(G2Sat));
    B(BSat) =   max( ( MinSatBal(BSat) .*(1-MinSat(BSat))  +  B(BSat)*fact .* MinSat(BSat) )/DecodeOptions.ColourBalance(3) ,  B(BSat));
    clear RSat G1Sat G2Sat BSat MinSat MinSatBal
end

%--- White Balance ---
if(DecodeOptions.EarlyWhiteBalance)
    % White Balance on RAW data before the demosaicing
    LensletImage(RstY:2:end,RstX:2:end) = R * DecodeOptions.ColourBalance(1);
    LensletImage(G1stY:2:end,G1stX:2:end) = G1 * DecodeOptions.ColourBalance(2);
    LensletImage(G2stY:2:end,G2stX:2:end) = G2 * DecodeOptions.ColourBalance(2);
    LensletImage(BstY:2:end,BstX:2:end) = B * DecodeOptions.ColourBalance(3);
elseif(DecodeOptions.CorrectSaturated)
    %Otherwise, the White Balance will be performed later
    LensletImage(RstY:2:end,RstX:2:end) = R;
    LensletImage(G1stY:2:end,G1stX:2:end) = G1;
    LensletImage(G2stY:2:end,G2stX:2:end) = G2;
    LensletImage(BstY:2:end,BstX:2:end) = B;
end
clear G1 G2 R B


if(DecodeOptions.WeightedDemosaic || DecodeOptions.WeightedInterp)
% P. DAVID: Computation of the lenslet centers and the lenslet mask
    [Belonging, MLCenters, Dist] = MicrolensBelonging(size(LensletImage,1),size(LensletImage,2), LensletGridModel);
    %Weights = (max(Dist(:))-Dist)/(max(Dist(:))-min(Dist(:)));
    clear Dist
    %%{
    Weights = cast(WhiteImage.*double(intmax('uint16')), 'uint16');
    Weights = demosaic(Weights, DecodeOptions.DemosaicOrder);
    Weights = cast(Weights, DecodeOptions.Precision);
    Weights = Weights ./  double(intmax('uint16'));
    Weights = RGB2YCbCr(Weights,1); Weights = Weights(:,:,1);
    %}
%{
%Display white image with centers
    MLCenters2=round(MLCenters);
    MLCenters2(MLCenters2(:,1)<=0 | MLCenters2(:,2)<=0,:)=[];
    MLCenters2(MLCenters2(:,1)>size(LensletImage,2) | MLCenters2(:,2)>size(LensletImage,1),:)=[];
    Temp=Weights;
    Temp(sub2ind(size(Temp),MLCenters2(:,2),MLCenters2(:,1)))=1;
    MLCentersDisp=Temp;
    Temp=Weights;
    Temp(sub2ind(size(Temp),MLCenters2(:,2),MLCenters2(:,1)))=0;
    MLCentersDisp=cat(3,MLCentersDisp,Temp,Temp);
    figure,imshow(MLCentersDisp,[]);
    clear Temp
%Display map of microlens belongings
    %figure,imshow(mod(Belonging,64),'colormap',colormap('prism'));caxis([0 64]);
    figure,imshow(mod(Belonging,1000),'colormap',rand(1000,3));caxis([0 1000]);
%}
end

if( nargout < 2 )
    clear WhiteImage
end

% Hot pixel correction
if(DecodeOptions.HotPixelCorrect)
    LensletImage = LFHotPixelCorrection(LensletImage,HotPixels);
end


LensletImage(~isfinite(LensletImage)) = 1;
LensletImage = max(0, LensletImage);%remove negative values
if(DecodeOptions.noClip)
    maxLum=max(LensletImage(:));
else
    %Apply soft clipping that keeps details in the highlights.
    maxLum=2^(-DecodeOptions.ExposureBias);
    LensletImage = softClip(LensletImage/maxLum,7)*maxLum;
end



%---Demosaic---
if DecodeOptions.WeightedDemosaic
    % P. DAVID: White lenslet image guided demosaicing
    LensletImage = weighted_demosaic( LensletImage, Weights.^10, Belonging, DecodeOptions );
else
    % This uses Matlab's demosaic, which is "gradient compensated". This likely has implications near
    % the edges of lenslet images, where the contrast is due to vignetting / aperture shape, and is not
    % a desired part of the image
    LensletImage = cast(LensletImage.*double(intmax('uint16'))./maxLum, 'uint16');
    LensletImage = demosaic(LensletImage, DecodeOptions.DemosaicOrder);
    LensletImage = cast(LensletImage, DecodeOptions.Precision);
    LensletImage = LensletImage ./  double(intmax('uint16')) .*maxLum;
end

DecodeOptions.NColChans = 3;

if( nargout >= 2 )
    DecodeOptions.NWeightChans = size(WhiteImage,3);
else
    DecodeOptions.NWeightChans = 0;
end

if( nargout > 3 )
    DebayerLensletImage = LensletImage;
end

%---Tranform to an integer-spaced grid---
fprintf('\nAligning image to lenslet array...');
InputSpacing = [LensletGridModel.HSpacing, LensletGridModel.VSpacing];
NewLensletSpacing = ceil(InputSpacing);
% Force even so hex shift is a whole pixel multiple
NewLensletSpacing = ceil(NewLensletSpacing/2)*2;
XformScale = NewLensletSpacing ./ InputSpacing;  % Notice the resized image will not be square

NewOffset = [LensletGridModel.HOffset, LensletGridModel.VOffset] .* XformScale;
RoundedOffset = round(NewOffset);
XformTrans =  RoundedOffset-NewOffset;

NewLensletGridModel = struct('HSpacing',NewLensletSpacing(1), 'VSpacing',NewLensletSpacing(2), ...
    'HOffset',RoundedOffset(1), 'VOffset',RoundedOffset(2), 'Rot',0, ...
    'UMax', LensletGridModel.UMax, 'VMax', LensletGridModel.VMax, 'Orientation', LensletGridModel.Orientation, ...
    'FirstPosShiftRow', LensletGridModel.FirstPosShiftRow);

%---Fix image rotation and scale---
RRot = LFRotz( LensletGridModel.Rot ); %counterclockwise rotation of angle Rot.

RScale = eye(3);
RScale(1,1) = XformScale(1);
RScale(2,2) = XformScale(2);
DecodeOptions.OutputScale(1:2) = XformScale;
DecodeOptions.OutputScale(3:4) = [1,2/sqrt(3)];  % hex sampling

RTrans = eye(3);
RTrans(end,1:2) = XformTrans;

% The following rotation can rotate parts of the lenslet image out of frame.
% todo[optimization]: attempt to keep these regions, offer greater user-control of what's kept
FixAll = maketform('affine', RRot*RScale*RTrans);%note mlp:apply rotation first, then scale, then translation.
NewSize = size(LensletImage(:,:,1)) .* XformScale(2:-1:1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MLP : Added Barycentric interpolation as in KAIST toolbox %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(DecodeOptions.ResampMethod,'barycentric'))
    rotateRaw=false;%true => rotate the row before barycentric interpolation / false => directy perform barycentric intepolation using rotated centers.
    if(rotateRaw)
        LensletGridModel2 = NewLensletGridModel;
        LensletImage = imtransform( LensletImage, FixAll, 'YData',[1 NewSize(1)], 'XData',[1 NewSize(2)]);
        if( nargout >= 2 )
            WhiteImage = imtransform( WhiteImage, FixAll, 'YData',[1 NewSize(1)], 'XData',[1 NewSize(2)]);
        end
    else
        LensletGridModel2 = LensletGridModel;
    end
    
    if( nargout >= 4 )
        CorrectedLensletImage = LensletImage;
    end
    
    fprintf('\nInterpolating subaperture images from raw data...');
    
    MaxSpacing = max(NewLensletGridModel.HSpacing, NewLensletGridModel.VSpacing);  % Enforce square output in s,t
    STVec = floor((-(MaxSpacing-1)/2)):ceil(((MaxSpacing-1)/2));
    FixLensletOffsets = inv(RRot*RScale);
    
    USize = LensletGridModel2.UMax;
    VSize = LensletGridModel2.VMax;
    HOffset = LensletGridModel2.HOffset;
    VOffset = LensletGridModel2.VOffset;
    HSpacing = LensletGridModel2.HSpacing;
    VSpacing = LensletGridModel2.VSpacing;
    HexShiftStart = LensletGridModel2.FirstPosShiftRow;
    
    [uu, vv] = meshgrid( 0:USize-1 , 0:VSize-1 );
    centerX = HOffset + HSpacing * uu;
    centerX(HexShiftStart:2:end,:) = centerX(HexShiftStart:2:end,:) + HSpacing/2;
    centerY = VOffset + VSpacing * vv;
    
    if(~rotateRaw)
        RotatedCenters = [centerX(:),centerY(:)] * RRot(1:2,1:2)';
        centerX(:) = RotatedCenters(:,1);
        centerY(:) = RotatedCenters(:,2);
    end
    
    VSizeNew = floor(VSize*3*sqrt(3)/2);
    USizeNew = USize*3;
    LF = ones(length(STVec), length(STVec), VSizeNew, USizeNew, DecodeOptions.NColChans + DecodeOptions.NWeightChans, DecodeOptions.Precision);
    
    first = true;
    for col = 1:length(STVec)
        for row = 1:length(STVec)
            if(~rotateRaw)
                vecFix = FixLensletOffsets * [STVec(row); STVec(col); 1];
                rowOffset = vecFix(1)/vecFix(3);
                colOffset = vecFix(2)/vecFix(3);
            else
                rowOffset = STVec(row);
                colOffset = STVec(col);
            end
            hex = zeros(2*VSize,2*USize,DecodeOptions.NColChans);
            %img = zeros(floor(2*VSize*VSpacing/HSpacing), 2*USize,DecodeOptions.NColChans);
            img = zeros(VSizeNew,USizeNew,DecodeOptions.NColChans);
            for ch = 1:DecodeOptions.NColChans
                InterpSlice = interpn(squeeze(LensletImage(:,:,ch)), centerY+colOffset, centerX+rowOffset, 'cubic');%'spline'
                hex(1:4:end,3-HexShiftStart:2:end,ch) =  InterpSlice(1:2:end,:);
                hex(3:4:end,HexShiftStart:2:end,ch) =  InterpSlice(2:2:end,:);
                if(first)
                    [img(:,:,ch), InterpData] = baryCentric(hex(:,:,ch),VSize/2,USize, HexShiftStart);
                    first = false;
                else
                   img(:,:,ch) = baryCentric(hex(:,:,ch),VSize/2,USize, HexShiftStart, InterpData);
                end
            end
            LF(col, row,:,:,1:DecodeOptions.NColChans) = cast(img,DecodeOptions.Precision);
        end
    end
    
else
    
    if DecodeOptions.WeightedInterp
    % P. DAVID: White lenslet image guided rotation
        LensletImage = lenslettransform(LensletImage, Weights.^10, LensletGridModel.Rot, XformTrans, XformScale, [1 NewSize(2)], [1 NewSize(1)], Belonging, MLCenters);
    else
        LensletImage = imtransform( LensletImage, FixAll, 'YData',[1 NewSize(1)], 'XData',[1 NewSize(2)]);
    end

    if( nargout >= 2 )
        WhiteImage = imtransform( WhiteImage, FixAll, 'YData',[1 NewSize(1)], 'XData',[1 NewSize(2)]);
    end

    if( nargout >= 4 )
        CorrectedLensletImage = LensletImage;
    end
    
    LF = SliceXYImage( NewLensletGridModel, LensletImage, WhiteImage, DecodeOptions );

    %---Correct for hex grid and resize to square u,v pixels---
    LFSize = size(LF);
    HexAspect = 2/sqrt(3);

    switch( DecodeOptions.ResampMethod )
        case 'fast'
            fprintf('\nResampling (1D approximation) to square u,v pixels');
            NewUVec = 0:1/HexAspect:(size(LF,4)+1);  % overshoot then trim
            NewUVec = NewUVec(1:ceil(LFSize(4)*HexAspect));
            OrigUSize = size(LF,4);
            LFSize(4) = length(NewUVec);
            %---Allocate dest and copy orig LF into it (memory saving vs. keeping both separately)---
            LF2 = zeros(LFSize, DecodeOptions.Precision);
            LF2(:,:,:,1:OrigUSize,:) = LF;
            LF = LF2;
            clear LF2

            if( DecodeOptions.DoDehex )
                ShiftUVec = -0.5+NewUVec;
                fprintf(' and removing hex sampling...');
            else
                ShiftUVec = NewUVec;
                fprintf('...');
            end
            for ColChan = 1:size(LF,5)
                CurUVec = ShiftUVec;
                for RowIter = 1:2
                    RowIdx = mod(NewLensletGridModel.FirstPosShiftRow + RowIter, 2) + 1;
                    ShiftRows = squeeze(LF(:,:,RowIdx:2:end,1:OrigUSize, ColChan));
                    SliceSize = size(ShiftRows);
                    SliceSize(4) = length(NewUVec);
                    ShiftRows = reshape(ShiftRows, [size(ShiftRows,1)*size(ShiftRows,2)*size(ShiftRows,3), size(ShiftRows,4)]);
                    ShiftRows = interp1( (0:size(ShiftRows,2)-1)', ShiftRows', CurUVec' )';
                    ShiftRows(isnan(ShiftRows)) = 0;
                    LF(:,:,RowIdx:2:end,:,ColChan) = reshape(ShiftRows,SliceSize);
                    CurUVec = NewUVec;
                end
            end
            clear ShiftRows
            DecodeOptions.OutputScale(3) = DecodeOptions.OutputScale(3) * HexAspect;

        case 'triangulation'
            fprintf('\nResampling (triangulation) to square u,v pixels');
            OldVVec = (0:size(LF,3)-1);
            OldUVec = (0:size(LF,4)-1) * HexAspect;

            NewUVec = (0:ceil(LFSize(4)*HexAspect)-1);
            NewVVec = (0:LFSize(3)-1);
            LFSize(4) = length(NewUVec);
            LF2 = zeros(LFSize, DecodeOptions.Precision);

            [Oldvv,Olduu] = ndgrid(OldVVec,OldUVec);
            [Newvv,Newuu] = ndgrid(NewVVec,NewUVec);
            if( DecodeOptions.DoDehex )
                fprintf(' and removing hex sampling...');
                FirstShiftRow = NewLensletGridModel.FirstPosShiftRow;
                Olduu(FirstShiftRow:2:end,:) = Olduu(FirstShiftRow:2:end,:) + HexAspect/2;
            else
                fprintf('...');
            end

            DT = delaunayTriangulation( Olduu(:), Oldvv(:) );  % use DelaunayTri in older Matlab versions
            [ti,bc] = pointLocation(DT, Newuu(:), Newvv(:));
            ti(isnan(ti)) = 1;

            for ColChan = 1:size(LF,5) 
                fprintf('.');
                for tidx= 1:LFSize(1)
                    for sidx= 1:LFSize(2)
                        CurUVSlice = squeeze(LF(tidx,sidx,:,:,ColChan));
                        triVals = CurUVSlice(DT(ti,:));
                        CurUVSlice = dot(bc',triVals')';
                        CurUVSlice = reshape(CurUVSlice, [length(NewVVec),length(NewUVec)]);

                        CurUVSlice(isnan(CurUVSlice)) = 0;
                        LF2(tidx,sidx, :,:, ColChan) = CurUVSlice;
                    end
                end
            end
            LF = LF2;
            clear LF2
            DecodeOptions.OutputScale(3) = DecodeOptions.OutputScale(3) * HexAspect;

        otherwise
            fprintf('\nNo valid dehex / resampling selected\n');
    end

    %---Resize to square s,t pixels---
    % Assumes only a very slight resampling is required, resulting in an identically-sized output light field
    if DecodeOptions.DoSquareST
        fprintf('\nResizing to square s,t pixels using 1D linear interp...');

        ResizeScale = DecodeOptions.OutputScale(1)/DecodeOptions.OutputScale(2);
        ResizeDim1 = 1;
        ResizeDim2 = 2;
        if( ResizeScale < 1 )
            ResizeScale = 1/ResizeScale;
            ResizeDim1 = 2;
            ResizeDim2 = 1;
        end

        OrigSize = size(LF, ResizeDim1);
        OrigVec = floor((-(OrigSize-1)/2):((OrigSize-1)/2));
        NewVec = OrigVec ./ ResizeScale;

        OrigDims = [1:ResizeDim1-1, ResizeDim1+1:5];

        UBlkSize = 32;
        USize = size(LF,4);
        LF = permute(LF,[ResizeDim1, OrigDims]);
        for UStart = 1:UBlkSize:USize
            UStop = UStart + UBlkSize - 1;
            UStop = min(UStop, USize);
            LF(:,:,:,UStart:UStop,:) = interp1(OrigVec, LF(:,:,:,UStart:UStop,:), NewVec);
            fprintf('.');
        end
        LF = ipermute(LF,[ResizeDim1, OrigDims]);
        LF(isnan(LF)) = 0;

        DecodeOptions.OutputScale(ResizeDim2) = DecodeOptions.OutputScale(ResizeDim2) * ResizeScale;
    end

end %endif ResampMethod == 'barycentric'

%---Trim s,t---
LF = LF(2:end-1,2:end-1, :,:, :);

%---Slice out LFWeight if it was requested---
if nargout >= 2
    LFWeight = LF(:,:,:,:,DecodeOptions.NColChans+1:end);
    for WChan=1:DecodeOptions.NWeightChans, LFWeight(:,:,:,:,WChan) = LFWeight(:,:,:,:,WChan)./max(max(max(max(LFWeight(:,:,:,:,WChan)))));end
    LF = LF(:,:,:,:,1:DecodeOptions.NColChans);
end

end


%------------------------------------------------------------------------------------------------------
function LF = SliceXYImage( LensletGridModel, LensletImage, WhiteImage, DecodeOptions )
% todo[optimization]: The SliceIdx and ValidIdx variables could be precomputed

fprintf('\nSlicing lenslets into LF...');

USize = LensletGridModel.UMax;
VSize = LensletGridModel.VMax;
MaxSpacing = max(LensletGridModel.HSpacing, LensletGridModel.VSpacing);  % Enforce square output in s,t
SSize = MaxSpacing + 1; % force odd for centered middle pixel -- H,VSpacing are even, so +1 is odd
TSize = MaxSpacing + 1;

LF = zeros(TSize, SSize, VSize, USize, DecodeOptions.NColChans + DecodeOptions.NWeightChans, DecodeOptions.Precision);

TVec = cast(floor((-(TSize-1)/2):((TSize-1)/2)), 'int16');
SVec = cast(floor((-(SSize-1)/2):((SSize-1)/2)), 'int16');
VVec = cast(0:VSize-1, 'int16');
UBlkSize = 32;
for UStart = 0:UBlkSize:USize-1  % note zero-based indexing
    UStop = UStart + UBlkSize - 1;
    UStop = min(UStop, USize-1);  
    UVec = cast(UStart:UStop, 'int16');
    
    [tt,ss,vv,uu] = ndgrid( TVec, SVec, VVec, UVec );
    
    %---Build indices into 2D image---
    LFSliceIdxX = LensletGridModel.HOffset + uu.*LensletGridModel.HSpacing + ss;
    LFSliceIdxY = LensletGridModel.VOffset + vv.*LensletGridModel.VSpacing + tt;
    
    HexShiftStart = LensletGridModel.FirstPosShiftRow;
    LFSliceIdxX(:,:,HexShiftStart:2:end,:) = LFSliceIdxX(:,:,HexShiftStart:2:end,:) + LensletGridModel.HSpacing/2;
    
    %---Lenslet mask in s,t and clip at image edges---
    CurSTAspect = DecodeOptions.OutputScale(1)/DecodeOptions.OutputScale(2);
    R = sqrt((cast(tt,DecodeOptions.Precision)*CurSTAspect).^2 + cast(ss,DecodeOptions.Precision).^2);
    ValidIdx = find(R < LensletGridModel.HSpacing/2 & ...
        LFSliceIdxX >= 1 & LFSliceIdxY >= 1 & LFSliceIdxX <= size(LensletImage,2) & LFSliceIdxY <= size(LensletImage,1) );
    
    %--clip -- the interp'd values get ignored via ValidIdx--
    LFSliceIdxX = max(1, min(size(LensletImage,2), LFSliceIdxX ));
    LFSliceIdxY = max(1, min(size(LensletImage,1), LFSliceIdxY ));
    
    %---
    LFSliceIdx = sub2ind([size(LensletImage,1),size(LensletImage,2)], cast(LFSliceIdxY,'int32'), ...
        cast(LFSliceIdxX,'int32'), ones(size(LFSliceIdxX),'int32'));
    
    tt = tt - min(tt(:)) + 1;
    ss = ss - min(ss(:)) + 1;
    vv = vv - min(vv(:)) + 1;
    uu = uu - min(uu(:)) + 1 + UStart;
    LFOutSliceIdx = sub2ind(size(LF), cast(tt,'int32'), cast(ss,'int32'), ...
        cast(vv,'int32'),cast(uu,'int32'), ones(size(ss),'int32'));
    
    %---
    for ColChan = 1:DecodeOptions.NColChans
        LF(LFOutSliceIdx(ValidIdx) + numel(LF(:,:,:,:,1)).*(ColChan-1)) = ...
            LensletImage( LFSliceIdx(ValidIdx) + numel(LensletImage(:,:,1)).*(ColChan-1) );
    end
    if DecodeOptions.NWeightChans ~= 0
        for WChan = 1:DecodeOptions.NWeightChans
            LF(LFOutSliceIdx(ValidIdx) + numel(LF(:,:,:,:,1)).*(DecodeOptions.NColChans + WChan - 1)) = ...
                WhiteImage( LFSliceIdx(ValidIdx) + numel(LensletImage(:,:,1)).*(WChan-1) );
        end
    end
    fprintf('.');
end
end



%------------------------------------------------------------------------------------------------------
%Soft clipping function (high value of R -> hard clipping)
function O = softClip(I,R)
b = exp(R);
O = log((1+b)./(1+b*exp(-R*I)))./log(1+b);
end