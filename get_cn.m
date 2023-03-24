function [ feature_image ] = get_cn( image, fparam, gparam )
w2crs=fparam.w2crs;
scaleNum=size(image,4);
feature_image=zeros(size(image,1),size(image,2),10,scaleNum);
if isempty(w2crs)
    % load the RGB to color name matrix if not in input
    load('w2crs');
end

for i=1:scaleNum
    feature_image(:,:,:,scaleNum) = im2c(single(hyperspectral2RGB(image(:,:,:,i))), w2crs, -2);
end