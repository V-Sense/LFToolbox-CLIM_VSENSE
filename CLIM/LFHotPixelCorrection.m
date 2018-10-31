%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original code from KAIST toolbox written by Donghyeon Cho.
%
% modified by Mikael Le Pendu
%  - adapted inputs/outputs for the LFtoolbox pipeline.
%
% Name   : hotPixelCorrection  
% Input  : img           - input image
%          mask_black    - hot pixel mask (from black images).
%
% Output : result        - corrected image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = LFHotPixelCorrection(img,mask_black)
    
    RawSize = size(img);
    
    mask_black(1:2,:) = 0;
    mask_black(:,1:2) = 0;
    mask_black(end-1:end,:) = 0;
    mask_black(:,end-1:end) = 0;
    index = find(mask_black == 1);
    
    [Y X] = ind2sub(RawSize,index);
    
    img(index) = 0;
    
    cnt = 1;
    sum = 0;
    for col=-1:1:1
        for row=-1:1:1
            w = exp(-sqrt(col^2+row^2));
            sum = sum + exp(-sqrt(col^2+row^2));
            ind = sub2ind(RawSize,Y+row,X+col);
            img(index) = img(index) + img(ind).*w;    
            cnt = cnt+1;            
        end
    end
    
    img(index) = img(index)./sum;
    result = img;
    
end