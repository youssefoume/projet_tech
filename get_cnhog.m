function [ feature_image ] = get_cnhog( im, fparam, gparam )
%extract fhog features using piotrs toolbox. Currently takes no parameters
%except hog-cell-size
if ~isfield(fparam, 'nOrients')
    fparam.nOrients = 9;
end
% bandNum=16;
[im_height, im_width, num_im_chan, num_images] = size(im);
feature_image = zeros(floor(im_height/gparam.cell_size), floor(im_width/gparam.cell_size), fparam.nDim, num_images, 'single');
hog_image= zeros(floor(im_height/gparam.cell_size), floor(im_width/gparam.cell_size), (fparam.nDim), 'single');
for k = 1:num_images
    for j=1: num_im_chan
        f=fhog(single(im(:,:,j,k)), gparam.cell_size, fparam.nOrients);
        hog_image(:,:,(j-1)*(fparam.nDim/num_im_chan)+1:(fparam.nDim/num_im_chan)*j)=f(:,:,1:end-1);
    end
    %the last dimension is all 0 so we can discard it
    feature_image(:,:,:,k) = hog_image;
end
end