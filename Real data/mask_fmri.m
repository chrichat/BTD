function [ dat1] = mask_fmri( unimg )
%Create masked data, reshape based on the minimum rectangle surrounding the
% brain in each slice. Keep the maskm for uncomplete tensor.
[x,y,z] = ind2sub(size(unimg(:,:,:,4)),find(unimg(:,:,:,4)>0));
dat1=unimg(min(x):max(x),min(y):max(y),min(z):max(z),:);
end

