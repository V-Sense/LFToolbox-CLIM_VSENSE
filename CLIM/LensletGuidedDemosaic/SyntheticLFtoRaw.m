function [Raw, LensletGridModel, Mask, LF_out] = SyntheticLFtoRaw(LF, theta)

    LF = double(LF);
    HexShiftStart = 1;
    %Make views hexagonal
    VVec = 1:sqrt(3)/2:size(LF,3);
    OldSize = size(LF);
    NewSize = OldSize;
    NewSize(3) = length(VVec);
    LF_out = zeros(NewSize);

    for ColChan = 1:size(LF_out,5)
        ScaleCols = squeeze(LF(:,:,:,:, ColChan));
        ScaleCols = permute(ScaleCols,[1, 2, 4, 3]);
        ScaleCols = reshape(ScaleCols, [size(ScaleCols,1)*size(ScaleCols,2)*size(ScaleCols,3), size(ScaleCols,4)]);
        ScaleCols = interp1( (1:size(ScaleCols,2))', ScaleCols', VVec' )';
        ScaleCols(isnan(ScaleCols)) = 0;
        ScaleCols = reshape(ScaleCols,[NewSize(1),NewSize(2),NewSize(4),NewSize(3)]);
        ScaleCols = permute(ScaleCols,[1, 2, 4, 3]);
        LF_out(:,:,:,:,ColChan) = ScaleCols;

        ShiftRows = squeeze(LF_out(:,:,HexShiftStart:2:end,:, ColChan));
        SliceSize = size(ShiftRows);
        ShiftRows = reshape(ShiftRows, [size(ShiftRows,1)*size(ShiftRows,2)*size(ShiftRows,3), size(ShiftRows,4)]);
        ShiftRows = interp1( (0:size(ShiftRows,2)-1)', ShiftRows', (0.5:size(ShiftRows,2)-0.5)' )';
        ShiftRows(isnan(ShiftRows)) = 0;
        LF_out(:,:,HexShiftStart:2:end,:,ColChan) = reshape(ShiftRows,SliceSize);
    end

    LF = LF_out;
    clear ShiftRows ScaleCols SliceSize

    %Make Raw
    SSize = size(LF,2);
    TSize = size(LF,1);
    USize = size(LF,4);
    VSize = size(LF,3);

    HSpacing = 2*ceil(max(SSize, TSize)/2);
    VSpacing = HSpacing*sqrt(3)/2;
    HOffset = (SSize+1)/2;
    VOffset = (TSize+1)/2;

    LF = single(LF);
    W = ceil(USize*HSpacing);
    H = ceil(VSize*VSpacing);
    Rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    corners = Rot*[1, 1, W, W, HOffset; 1, H, 1, H, VOffset];
    W = ceil(max(corners(1,1:end-1)) - min(corners(1,1:end-1)));
    H = ceil(max(corners(2,1:end-1)) - min(corners(2,1:end-1)));
    
    LensletGridModel.HSpacing = HSpacing;
    LensletGridModel.VSpacing = VSpacing;
    LensletGridModel.HOffset =  corners(1,end) - min(corners(1,1:end-1)) + 1;
    LensletGridModel.VOffset =  corners(2,end) - min(corners(2,1:end-1)) + 1;
    LensletGridModel.Rot = theta;
    LensletGridModel.Orientation = 'horz';
    LensletGridModel.FirstPosShiftRow = HexShiftStart;
    LensletGridModel.UMax = USize;
    LensletGridModel.VMax = VSize;
    
    TVec = floor((-(TSize-1)/2):((TSize-1)/2));
    SVec = floor((-(SSize-1)/2):((SSize-1)/2));
    VVec = 0:VSize-1;
    UBlkSize = 32;
    Raw = zeros(H, W, size(LF,5),'single');
    Mask = zeros(H,W,'single');

    for UStart = 0:UBlkSize:USize-1  % note zero-based indexing
        UStop = UStart + UBlkSize - 1;
        UStop = min(UStop, USize-1);  
        UVec = UStart:UStop;

        [tt,ss,vv,uu] = ndgrid( TVec, SVec, VVec, UVec );

        %---Build indices into 2D image---
        LFSliceIdxX = uu.*LensletGridModel.HSpacing + ss;
        LFSliceIdxY = vv.*LensletGridModel.VSpacing + tt;

        LFSliceIdxX(:,:,LensletGridModel.FirstPosShiftRow:2:end,:) = LFSliceIdxX(:,:,LensletGridModel.FirstPosShiftRow:2:end,:) + HOffset;
        
        LFSliceIdx = Rot*[LFSliceIdxX(:)';LFSliceIdxY(:)'];
        LFSliceIdxX = round(LensletGridModel.HOffset + reshape(LFSliceIdx(1,:),size(uu)));
        LFSliceIdxY = round(LensletGridModel.VOffset + reshape(LFSliceIdx(2,:),size(uu)));
        clear LFSliceIdx
        
        %---Lenslet mask in s,t and clip at image edges---
        R = sqrt(tt.^2 + ss.^2);
        ValidIdx = find(R < (SSize+1)/2 & ...
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
        
        VignetOutSliceIdx = sub2ind(size(R), cast(tt,'int32'), cast(ss,'int32'), ...
            cast(vv,'int32'),cast(uu - UStart,'int32'), ones(size(ss),'int32'));
        %---          
        vignet = cos(R(VignetOutSliceIdx(ValidIdx))*0.2).^4;
        for ColChan = 1:size(LF,5) 
            Raw( LFSliceIdx(ValidIdx) + numel(Raw(:,:,1)).*(ColChan-1) ) = ...
                LF(LFOutSliceIdx(ValidIdx) + numel(LF(:,:,:,:,1)).*(ColChan-1)).*vignet;
        end
        
        Mask(LFSliceIdx(ValidIdx)) = 255.*vignet;
        fprintf('.');
    end
    
%     for ColChan = 1:size(Raw,3)
%             Raw(:,:,ColChan) = awgn(Raw(:,:,ColChan),20,'measured');
%     end
    Raw = min(255, max(Raw,0));
end