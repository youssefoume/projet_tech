% This function implements the BACF tracker.

function [results] = MHT_optimized(params)

%   Setting parameters for local use.
search_area_scale   = params.search_area_scale;
output_sigma_factor = params.output_sigma_factor;
learning_rate       = params.learning_rate;
filter_max_area     = params.filter_max_area;
nScales             = params.number_of_scales;
scale_step          = params.scale_step;
interpolate_response = params.interpolate_response;

features    = params.t_features;
video_path  = params.video_path;
s_frames    = params.s_frames;
pos         = floor(params.init_pos);
target_sz   = floor(params.wsize);

visualization  = params.visualization;
num_frames     = params.no_fram;
init_target_sz = target_sz;

%set the feature ratio to the feature-cell size
featureRatio = params.t_global.cell_size;
search_area = prod(init_target_sz / featureRatio * search_area_scale);
feature_selection=params.feature_selection;
% when the number of cells are small, choose a smaller cell size
if isfield(params.t_global, 'cell_selection_thresh')
    if search_area < params.t_global.cell_selection_thresh * filter_max_area
        %         params.t_global.cell_size = min(featureRatio, max(1, ceil(sqrt(prod(init_target_sz * search_area_scale)/(params.t_global.cell_selection_thresh * filter_max_area)))));
        
        featureRatio = params.t_global.cell_size;
        search_area = prod(init_target_sz / featureRatio * search_area_scale);
    end
end

global_feat_params = params.t_global;

if search_area > filter_max_area
    currentScaleFactor = sqrt(search_area / filter_max_area);
else
    currentScaleFactor = 1.0;
end

% target size at the initial scale
base_target_sz = target_sz / currentScaleFactor;
% window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        sz = floor( base_target_sz * search_area_scale);     % proportional area, same aspect ratio as the target
    case 'square'
        sz = repmat(sqrt(prod(base_target_sz * search_area_scale)), 1, 2); % square area, ignores the target aspect ratio
    case 'fix_padding'
        sz = base_target_sz + sqrt(prod(base_target_sz * search_area_scale) + (base_target_sz(1) - base_target_sz(2))/4) - sum(base_target_sz)/2; % const padding
    otherwise
        error('Unknown "params.search_area_shape". Must be ''proportional'', ''square'' or ''fix_padding''');
end

% set the size to exactly match the cell size
sz = round(sz / featureRatio) * featureRatio;
use_sz = floor(sz/featureRatio);

% construct the label function- correlation output, 2D gaussian function,
% with a peak located upon the target
output_sigma = sqrt(prod(floor(base_target_sz/featureRatio))) * output_sigma_factor;
rg           = circshift(-floor((use_sz(1)-1)/2):ceil((use_sz(1)-1)/2), [0 -floor((use_sz(1)-1)/2)]);
cg           = circshift(-floor((use_sz(2)-1)/2):ceil((use_sz(2)-1)/2), [0 -floor((use_sz(2)-1)/2)]);
[rs, cs]     = ndgrid( rg,cg);
y            = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
yf           = fft2(y); %   FFT of y.

if interpolate_response == 1
    interp_sz = use_sz * featureRatio;
else
    interp_sz = use_sz;
end

% construct cosine window
cos_window = single(hann(use_sz(1))*hann(use_sz(2))');

% Calculate feature dimension
try
    im=imread([video_path '/HSI/' s_frames{1}]);
catch
    try
        im = imread(s_frames{1});
    catch
        %disp([video_path '/' s_frames{1}])
        im = imread([video_path '/' s_frames{1}]);
    end
end
% im=im(33:32+1024,1:2048);
im = Demosaicing_img(im);
%im=X2Cube(im);
im=im./max(im(:));
if size(im,3) == 3
    if all(all(im(:,:,1) == im(:,:,2)))
        colorImage = false;
    else
        colorImage = true;
    end
else
    colorImage = false;
end

% compute feature dimensionality
feature_dim = 0;
im = imread([video_path '/HSI/' s_frames{1}]);
% im=im(33:32+1024,1:2048);
im = Demosaicing_img(im);
%im=X2Cube(im);
im=im./max(im(:));
fIndex=[];
for n = 1:length(features)
    
    if strcmp(features{n}.feature,'material')
        features{n}.fparams.useForMaterial = true;
        % extract training sample image region
        pixels = get_pixels(im,pos,init_target_sz,init_target_sz);
        %         load spectralLibNew.mat;
        load(features{n}.fparams.endFileName);
        lib=max(lib,eps);
        lib=selectMaterial(pixels,lib);
        feature_dim=feature_dim+size(lib,1);
        features{n}.fparams.lib=lib;
        features{n}.fparams.nDim=size(lib,1);
        fIndex=[fIndex,size(lib,1)];
        
    end
    
    if strcmp(features{n}.feature,'hog')
        features{n}.fparams.useForHog = true;
        feature_dim=feature_dim+ features{n}.fparams.nDim;
        fIndex=[fIndex,features{n}.fparams.nDim];
    end
    if strcmp(features{n}.feature,'spectral')
        features{n}.fparams.useForSpectral = true;
        feature_dim=feature_dim+ size(im,3);
        fIndex=[fIndex,size(im,3)];
        
    end
    if strcmp(features{n}.feature,'cn')
        features{n}.fparams.useForCN = true;
        feature_dim=feature_dim+ features{n}.fparams.nDim;
        load(features{n}.fparams.w2c)
        features{n}.fparams.w2crs=w2crs;
        fIndex=[fIndex,features{n}.fparams.nDim];
    end
end
if feature_selection
    fw=ones(1,length(fIndex))./length(fIndex);
    fIndex=[0 cumsum(fIndex)];
    
end

% if size(im,3) > 1 && colorImage == false
%     im = im(:,:,1);
% end

if nScales > 0
    scale_exp = (-floor((nScales-1)/2):ceil((nScales-1)/2));
    scaleFactors = scale_step .^ scale_exp;
    min_scale_factor = scale_step ^ ceil(log(max(5 ./ sz)) / log(scale_step));
    max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));
end

if interpolate_response >= 3
    % Pre-computes the grid that is used for socre optimization
    ky = circshift(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), [1, -floor((use_sz(1) - 1)/2)]);
    kx = circshift(-floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2), [1, -floor((use_sz(2) - 1)/2)])';
    newton_iterations = params.newton_iterations;
end

% initialize the projection matrix (x,y,h,w)
rect_position = zeros(num_frames, 4);
time = 0;

% allocate memory for multi-scale tracking
multires_pixel_template = zeros(sz(1), sz(2), size(im,3), nScales);
small_filter_sz = floor(base_target_sz/featureRatio);

loop_frame = 1;
for frame = 1:numel(s_frames)
    %load image
    try
        im = imread([video_path '/HSI/' s_frames{frame}]);
    catch
        try
            im = imread([s_frames{frame}]);
        catch
            im = imread([video_path '/' s_frames{frame}]);
        end
    end
    disp(['Running frame ',num2str(frame)]);
    im = Demosaicing_img(im);
%     if frame==10
%         break;
%     end
   % im=X2Cube(im);
    im=im./max(im(:));
    %     if size(im,3) > 1 && colorImage == false
    %         im = im(:,:,1);
    %     end
    
    tic();
    
    %do not estimate translation and scaling on the first frame, since we
    %just want to initialize the tracker there
    if frame > 1
        for scale_ind = 1:nScales
            multires_pixel_template(:,:,:,scale_ind) = ...
                get_pixels(im, pos, round(sz*currentScaleFactor*scaleFactors(scale_ind)), sz);
        end
        xtf = fft2(bsxfun(@times,get_features(multires_pixel_template,features,global_feat_params),cos_window));
        if feature_selection
            responset = bsxfun(@times, conj(g_f), xtf);
            for t=1:length(fIndex)-1
                responset(:,:,fIndex(t)+1:fIndex(t+1),:)=responset(:,:,fIndex(t)+1:fIndex(t+1),:)*fw(t);
            end
            %             responset(:,:,1:size(lib,1),:)=responset(:,:,1:size(lib,1),:)*fw(1);
            %             responset(:,:,size(lib,1)+1:end,:)=responset(:,:,size(lib,1)+1:end,:)*fw(2);
            responsef=responset;
            responsef = permute(sum(responsef, 3), [1 2 4 3]);
        else
            responsef = permute(sum(bsxfun(@times, conj(g_f), xtf), 3), [1 2 4 3]);
        end
        % if we undersampled features, we want to interpolate the
        % response so it has the same size as the image patch
        if interpolate_response == 2
            % use dynamic interp size
            interp_sz = floor(size(y) * featureRatio * currentScaleFactor);
        end
        responsef_padded = resizeDFT2(responsef, interp_sz);
        
        % response in the spatial domain
        response = ifft2(responsef_padded, 'symmetric');
        
        % find maximum peak
        if interpolate_response == 3
            error('Invalid parameter value for interpolate_response');
        elseif interpolate_response == 4
            [disp_row, disp_col, sind,row,col] = resp_newton(response, responsef_padded, newton_iterations, ky, kx, use_sz);
        else
            [row, col, sind] = ind2sub(size(response), find(response == max(response(:)), 1));
            disp_row = mod(row - 1 + floor((interp_sz(1)-1)/2), interp_sz(1)) - floor((interp_sz(1)-1)/2);
            disp_col = mod(col - 1 + floor((interp_sz(2)-1)/2), interp_sz(2)) - floor((interp_sz(2)-1)/2);
        end
        % calculate translation
        switch interpolate_response
            case 0
                translation_vec = round([disp_row, disp_col] * featureRatio * currentScaleFactor * scaleFactors(sind));
            case 1
                translation_vec = round([disp_row, disp_col] * currentScaleFactor * scaleFactors(sind));
            case 2
                translation_vec = round([disp_row, disp_col] * scaleFactors(sind));
            case 3
                translation_vec = round([disp_row, disp_col] * featureRatio * currentScaleFactor * scaleFactors(sind));
            case 4
                translation_vec = round([disp_row, disp_col] * featureRatio * currentScaleFactor * scaleFactors(sind));
        end
        % set the scale
        temp_scale=currentScaleFactor;
        currentScaleFactor = currentScaleFactor * scaleFactors(sind);
        % adjust to make sure we are not to large or to small
        if currentScaleFactor < min_scale_factor
            currentScaleFactor = min_scale_factor;
        elseif currentScaleFactor > max_scale_factor
            currentScaleFactor = max_scale_factor;
        end
        
        % update position
        old_pos = pos;
        pos = pos + translation_vec;
        target_sz = floor(base_target_sz * currentScaleFactor);
        rect_pre=[pos([2,1]) - floor(target_sz([2,1])/2), target_sz([2,1])];
        if feature_selection
            disp_row_=zeros(1,length(fIndex)-1);
            disp_col_=disp_row_;
            sind_=disp_row_;
            row_=disp_row_;
            col_=disp_row_;
            responsef_=zeros(size(responset,1),size(responset,2),size(responset,4));
            responsef_padded_=responsef_;
            response_=responsef_;
            translation_vec_=zeros(2, length(fIndex)-1);
            currentScaleFactor_=zeros(length(fIndex)-1, 1);
            rect_pre_=zeros(4, length(fIndex)-1);
            overlaps_=zeros(length(fIndex)-1, 1);
            distances_=overlaps_;
            fw_overlap=overlaps_;
            fw_distance=overlaps_;
            fw_self=overlaps_;
            fw_=overlaps_;
        
            for t=1:length(fIndex)-1
                responsef_(:,:,:,t)=permute(sum(responset(:,:,fIndex(t)+1:fIndex(t+1),:), 3), [1 2 4 3]);
                responsef_padded_(:,:,:,t) = resizeDFT2(responsef_(:,:,:,t), interp_sz);
                response_(:,:,:,t)  = ifft2(responsef_padded_(:,:,:,t) , 'symmetric');
                if interpolate_response == 3
                    error('Invalid parameter value for interpolate_response');
                elseif interpolate_response == 4
                    [disp_row_(t), disp_col_(t), sind_(t),row_(t),col_(t)] = resp_newton(response_(:,:,:,t), responsef_padded_(:,:,:,t), newton_iterations, ky, kx, use_sz);
                else
                    [row_(t), col_(t), sind_(t)] = ind2sub(size(response_(:,:,:,t)), find(response_(:,:,:,t)== max(reshape(response_(:,:,:,t),[],1), 1)));
                    disp_row_(t) = mod(row_(t) - 1 + floor((interp_sz(1)-1)/2), interp_sz(1)) - floor((interp_sz(1)-1)/2);
                    disp_col_(t) = mod(col_(t) - 1 + floor((interp_sz(2)-1)/2), interp_sz(2)) - floor((interp_sz(2)-1)/2);
                end
                
                
                % calculate translation
                switch interpolate_response
                    case 0
                        translation_vec_(:,t)= round([disp_row_(t), disp_col_(t)] * featureRatio * temp_scale * scaleFactors(sind_(t)));
                    case 1
                        translation_vec_(:,t) = round([disp_row_(t), disp_col_(t)] * temp_scale * scaleFactors(sind_(t)));
                    case 2
                        translation_vec_(:,t) = round([disp_row_(t), disp_col_(t)] * scaleFactors(sind_m));
                    case 3
                        translation_vec_(:,t) = round([disp_row_(t), disp_col_(t)] * featureRatio * temp_scale * scaleFactors(sind_(t)));
                    case 4
                        translation_vec_(:,t) = round([disp_row_(t), disp_col_(t)] * featureRatio * temp_scale * scaleFactors(sind_(t)));
                end
                currentScaleFactor_(t) = temp_scale * scaleFactors(sind_(t));
                % adjust to make sure we are not to large or to small
                if currentScaleFactor_(t) < min_scale_factor
                    currentScaleFactor_(t) = min_scale_factor;
                elseif currentScaleFactor_(t)> max_scale_factor
                    currentScaleFactor_(t) = max_scale_factor;
                end
                pos_=old_pos+ translation_vec_(:,t)';
                target_sz_= floor(base_target_sz * currentScaleFactor_(t));
                rect_pre_(:,t)=[pos_([2,1]) - floor(target_sz_([2,1])/2), target_sz_([2,1])];
                [overlaps_(t),distances_(t)] =compute_relaibitlity(rect_pre_(:,t)', rect_pre);
                fw_overlap(t)=[exp(-(1-overlaps_(t))^2)];
                fw_distance(t)=[exp(-(distances_(t))^2)];
                fw_self(t)=[exp(-sum(translation_vec_(:,t).^2)/(2*mean(target_sz_)))];
                fw_(t)=fw_overlap(t)+fw_distance(t)+fw_self(t);
            end
            
            fw_t=fw_'./sum(fw_);
            fw=(1-learning_rate)*fw+learning_rate*fw_t;
                      
        end
        
        
    end
    
    % extract training sample image region
    pixels = get_pixels(im,pos,round(sz*currentScaleFactor),sz);
    
    % extract features and do windowing
    xf = fft2(bsxfun(@times,get_features(pixels,features,global_feat_params),cos_window));
    
    if (frame == 1)
        model_xf = xf;
    else
        model_xf = ((1 - learning_rate) * model_xf) + (learning_rate * xf);
    end
    
    g_f = single(zeros(size(xf)));
    h_f = g_f;
    l_f = g_f;
    mu    = 1;
    betha = 10;
    mumax = 10000;
    i = 1;
    temp_xf=model_xf;
    if (feature_selection)
        model_xf(:,:,1:size(lib,1))=temp_xf(:,:,1:size(lib,1))*fw(1);
        model_xf(:,:,size(lib,1)+1:end)=temp_xf(:,:,size(lib,1)+1:end)*fw(2);
    end
    T = prod(use_sz);
    S_xx = sum(conj(model_xf) .* model_xf, 3);
    params.admm_iterations = 2;
    %   ADMM
    while (i <= params.admm_iterations)
        %   solve for G- please refer to the paper for more details
        B = S_xx + (T * mu);
        S_lx = sum(conj(model_xf) .* l_f, 3);
        S_hx = sum(conj(model_xf) .* h_f, 3);
        g_f = (((1/(T*mu)) * bsxfun(@times, yf, model_xf)) - ((1/mu) * l_f) + h_f) - ...
            bsxfun(@rdivide,(((1/(T*mu)) * bsxfun(@times, model_xf, (S_xx .* yf))) - ((1/mu) * bsxfun(@times, model_xf, S_lx)) + (bsxfun(@times, model_xf, S_hx))), B);
        
        
        
        %   solve for H
        h = (T/((mu*T)+ params.admm_lambda))* ifft2((mu*g_f) + l_f);
        [sx,sy,h] = get_subwindow_no_window(h, floor(use_sz/2) , small_filter_sz);
        t = single(zeros(use_sz(1), use_sz(2), size(h,3)));
        t(sx,sy,:) = h;

        h_f = fft2(t);
        
        %   update L
        l_f = l_f + (mu * (g_f - h_f));
        
        %   update mu- betha = 10.
        mu = min(betha * mu, mumax);
        i = i+1;
    end
    model_xf=temp_xf;
    
    target_sz = floor(base_target_sz * currentScaleFactor);
    
    %save position and calculate FPS
    rect_position(loop_frame,:) = [pos([2,1]) - floor(target_sz([2,1])/2), target_sz([2,1])];
    
    time = time + toc();
    
    %visualization
    if visualization == 1
        rect_position_vis = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        im_to_show = double(hyper2im(im))/255;
        if size(im_to_show,3) == 1
            im_to_show = repmat(im_to_show, [1 1 3]);
        end
        if frame == 1
            fig_handle = figure('Name', 'Tracking');
            imagesc(im_to_show);
            hold on;
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
            text(10, 10, int2str(frame), 'color', [0 1 1]);
            hold off;
            axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
        else
            resp_sz = round(sz*currentScaleFactor*scaleFactors(scale_ind));
            xs = floor(old_pos(2)) + (1:resp_sz(2)) - floor(resp_sz(2)/2);
            ys = floor(old_pos(1)) + (1:resp_sz(1)) - floor(resp_sz(1)/2);
            sc_ind = floor((nScales - 1)/2) + 1;
            
            figure(fig_handle);
            imagesc(im_to_show);
            hold on;
            resp_handle = imagesc(xs, ys, fftshift(response(:,:,sc_ind))); colormap hsv;
            alpha(resp_handle, 0.2);
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
            text(20, 30, ['# Frame : ' int2str(loop_frame) ' / ' int2str(num_frames)], 'color', [1 0 0], 'BackgroundColor', [1 1 1], 'fontsize', 16);
            text(20, 60, ['FPS : ' num2str(1/(time/loop_frame))], 'color', [1 0 0], 'BackgroundColor', [1 1 1], 'fontsize', 16);
            
            hold off;
        end
        drawnow
    end
    loop_frame = loop_frame + 1;
end
%   save resutls.
fps = loop_frame / time;
results.type = 'rect';
results.res = rect_position;
results.fps = fps;
