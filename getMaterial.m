
function [X_hat_tv_i,param]=getMaterial(image,param,~)
endLib=param.lib;
endNum=param.nDim;
scaleNum=size(image,4);
X_hat_tv_i=zeros(size(image,1),size(image,2),endNum,scaleNum);
for i=1:scaleNum
    x_n=reshape(image(:,:,:,i),size(image,1)*size(image,2),size(image,3));
%     tmp=hyperFcls(x_n',endLib');
    [~, tmp] = simplex_project(x_n',endLib');
%         [~, tmp]=ASU(x_n',endLib');
%     [A_est, S_est, time] = HyperCSI(X,N)
    X_hat_tv_i(:,:,:,i)=reshape(tmp',size(image,1),size(image,2),endNum);
% constraint = 'sto';
% hesscomp = 'blockwise';
% [tmp,Cstd,Mem] = ipls(x_n',endLib', 'hesscomp', hesscomp, 'constraint', constraint);
%  X_hat_tv_i(:,:,:,i)=reshape(tmp',size(image,1),size(image,2),endNum);

%     tmp=x_n*pinv(endLib);
% %     [tmp,res,rmse_i] = sunsal_tv(endLib',x_n','MU',0.05,'POSITIVITY','yes','ADDONE','yes', ...
% %         'LAMBDA_1',lambda,'LAMBDA_TV', lambda_TV, 'TV_TYPE','iso',...
% %         'IM_SIZE',[size(image,1),size(image,2)],'AL_ITERS',100,  'VERBOSE','yes');
%     X_hat_tv_i(:,:,:,i)=reshape(tmp,size(image,1),size(image,2),endNum);
end
 
end