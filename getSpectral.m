function [X_hat_tv_i,param]=getSpectral(image,param,~)
ndim=param.nDim;
scaleNum=size(image,4);
X_hat_tv_i=zeros(size(image,1),size(image,2),ndim,scaleNum);
for i=1:scaleNum
    sz=size(image(:,:,:,i));
    X_hat_tv_i(:,:,:,i)=image(:,:,:,i)-reshape(repmat(mean(reshape(image(:,:,:,i),[],sz(3)),1),size(image,1)*size(image,2),1),sz);
end
end