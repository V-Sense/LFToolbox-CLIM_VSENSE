function [Raw, WhiteImage] = InvSliceXYImage(M,N,LensletGridModel, LF, DecodeOptions )
fprintf('\nMapping LF into a lenslet image...');

USize = LensletGridModel.UMax;
VSize = LensletGridModel.VMax;
MaxSpacing = max(LensletGridModel.HSpacing, LensletGridModel.VSpacing);
SSize = MaxSpacing + 1;
TSize = MaxSpacing + 1;

Raw = -ones(M,N, DecodeOptions.NColChans,DecodeOptions.Precision);
WhiteImage = zeros(M,N, DecodeOptions.NWeightChans,DecodeOptions.Precision);

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
        LFSliceIdxX >= 1 & LFSliceIdxY >= 1 & LFSliceIdxX <= size(Raw,2) & LFSliceIdxY <= size(Raw,1) );
    
    %--clip -- the interp'd values get ignored via ValidIdx--
    LFSliceIdxX = max(1, min(size(Raw,2), LFSliceIdxX ));
    LFSliceIdxY = max(1, min(size(Raw,1), LFSliceIdxY ));
    
    %---
    LFSliceIdx = sub2ind([size(Raw,1), size(Raw,2)], cast(LFSliceIdxY,'int32'), ...
        cast(LFSliceIdxX,'int32'), ones(size(LFSliceIdxX),'int32'));
    
    tt = tt - min(tt(:)) + 1;
    ss = ss - min(ss(:)) + 1;
    vv = vv - min(vv(:)) + 1;
    uu = uu - min(uu(:)) + 1 + UStart;
    LFOutSliceIdx = sub2ind(size(LF), cast(tt,'int32'), cast(ss,'int32'), ...
        cast(vv,'int32'),cast(uu,'int32'), ones(size(ss),'int32'));
    
    %---
    for ColChan = 1:DecodeOptions.NColChans 
        Raw( LFSliceIdx(ValidIdx) + numel(Raw(:,:,1)).*(ColChan-1) ) = ...
            LF(LFOutSliceIdx(ValidIdx) + numel(LF(:,:,:,:,1)).*(ColChan-1));
    end
    if DecodeOptions.NWeightChans ~= 0
        WhiteImage( LFSliceIdx(ValidIdx) ) = ...
            LF(LFOutSliceIdx(ValidIdx) + numel(LF(:,:,:,:,1)).*(DecodeOptions.NColChans));
    end
    fprintf('.');
end

        mask = (Raw(:,:,1) == -1) | (Raw(:,:,2) == -1) | (Raw(:,:,3) == -1);
        Raw(:,:,1) = regionfill(Raw(:,:,1), mask);
        Raw(:,:,2) = regionfill(Raw(:,:,2), mask);
        Raw(:,:,3) = regionfill(Raw(:,:,3), mask);
end