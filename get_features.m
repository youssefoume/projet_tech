
function [feature_pixels, support_sz,features] = get_features(image, features, gparams, fg_size )

if ~ iscell(features)
    features = {features};
end;

[im_height, im_width, num_im_chan, num_images] = size(image);

colorImage = num_im_chan == 3;


%compute total dimension of all features
tot_feature_dim = 0;
for n = 1:length(features)
    
   if strcmp(features{n}.feature,'material')
        features{n}.fparams.useForMaterial = true;
    end
    
    if strcmp(features{n}.feature,'hog')
        features{n}.fparams.useForHog = true;
    end
    if strcmp(features{n}.feature,'spectral')
        features{n}.fparams.useForSpectral = true;
    end
    tot_feature_dim = tot_feature_dim + features{n}.fparams.nDim;

    
end;

if nargin < 4 || isempty(fg_size)
    if gparams.cell_size == -1
        fg_size = size(features{1}.getFeature(image,features{1}.fparams,gparams));
    else
        fg_size = [floor(im_height/gparams.cell_size), floor(im_width/gparams.cell_size)];
    end
end

% temporary hack for fixing deep features
if gparams.cell_size == -1
    cf = features{1};
    if (cf.fparams.useForColor && colorImage) || (cf.fparams.useForGray && ~colorImage)
        [feature_pixels, support_sz] = cf.getFeature(image,cf.fparams,gparams);
    end;
else
    %compute the feature set
    feature_pixels = zeros(fg_size(1),fg_size(2),tot_feature_dim, num_images, 'single');
    
    currDim = 1;
    for n = 1:length(features)
        cf = features{n};
        if  (isfield(features{n}.fparams,'useForMaterial')&&features{n}.fparams.useForMaterial)||...
                isfield(features{n}.fparams,'useForSpectral')&&features{n}.fparams.useForSpectral 
%             image=imresize(image,fg_size);
            [fea,param]=cf.getFeature(imresize(image,fg_size),cf.fparams,gparams);
            features{n}.fparams=param;
%             fea=imresize(fea,fg_size);
            feature_pixels(:,:,currDim:(currDim+cf.fparams.nDim-1),:) = fea;
            currDim = currDim + cf.fparams.nDim;
            continue
        end
        
        if  (isfield(features{n}.fparams,'useForCN'))
%             image=imresize(image,fg_size);
            [fea]=cf.getFeature(imresize(image,fg_size),cf.fparams,gparams);
%             fea=imresize(fea,fg_size);
            feature_pixels(:,:,currDim:(currDim+cf.fparams.nDim-1),:) = fea;
            currDim = currDim + cf.fparams.nDim;
            continue
        end
        
        feature_pixels(:,:,currDim:(currDim+cf.fparams.nDim-1),:) =cf.getFeature(image,cf.fparams,gparams);
        currDim = currDim + cf.fparams.nDim;
    end;
    support_sz = [im_height, im_width];
end

end