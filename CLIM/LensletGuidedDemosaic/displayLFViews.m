function [ Image, Video ] = displayLFViews( LF )
%DISPLAYLFVIEWS Summary of this function goes here
%   Detailed explanation goes here
LF = LF(:,:,:,:,1:3);
LF = cast(LF, 'double');
LF = (LF - min(LF(:)))/(max(LF(:)) - min(LF(:))); 
LF = uint8(255*LF);

LF1 = permute(LF, [3 1 4 2 5]);
Image = reshape(LF1, size(LF1,1)*size(LF1,2),size(LF1,3)*size(LF1,4),size(LF1,5));

figure;imshow(Image);
Video = zeros(size(LF,3), size(LF,4),3, size(LF,1)*size(LF,2), 'uint8');

i=0;
for u = 1:size(LF,1)
    for v = 1:size(LF,2)
        i = i+1;
        Video(:,:,:,i) = squeeze(LF(u,v,:,:,:));
    end
end

%% CIRCLE
% N_frames = 100;
% pitch = 25;
% r = floor(min(size(LF,1), size(LF,2))/2)-1;
% center_u = ceil(size(LF,1)/2);
% center_v = ceil(size(LF,2)/2);
% Video = zeros(size(LF,3), size(LF,4),3, N_frames, 'uint8');
% 
% for i=0:N_frames-1
%     u = center_u + round(r*cos(2*pi*i/pitch));
%     v = center_v + round(r*sin(2*pi*i/pitch));
%     Video(:,:,:,i+1) = squeeze(LF(u,v,:,:,:));
% end

%% RECTANGLE
% N_laps = 4;
% N_framesperlap = (size(LF,1) + size(LF,2) - 6)*2;
% Video = zeros(size(LF,3), size(LF,4),3, N_framesperlap, 'uint8');
% 
% k = 0;
% for i=1:size(LF,1)-2
%     k = k+1;
%     Video(:,:,:,k) = squeeze(LF(1+i,2,:,:,:));
% end
% 
% for i=1:size(LF,2)-4
%     k = k+1;
%     Video(:,:,:,k) = squeeze(LF(size(LF,1)-1,2+i,:,:,:));
% end
% 
% for i=1:size(LF,1)-2
%     k = k+1;
%     Video(:,:,:,k) = squeeze(LF(size(LF,1)-i,size(LF,2)-1,:,:,:));
% end
% 
% for i=1:size(LF,2)-4
%     k = k+1;
%     Video(:,:,:,k) = squeeze(LF(2,size(LF,2)-1-i,:,:,:));
% end
% 
% Video = repmat(Video,1,1,1,N_laps);

implay(Video,30);
end

