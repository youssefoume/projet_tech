function [X_hat_tv_i, param]=getDymaticMaterial(image2,param,gparam)
%%%% image1£º previous frame
%%%image 2: current  frame 
%%% param contains the lib for each scale 
%
endNum=param.nDim;
scaleNum=size(image2,4);
image1=gparam.preFrame;
X_hat_tv_i=zeros(size(image2,1),size(image2,2),endNum,scaleNum);
%%%%%for dynamic  unmixing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lS = 1000000;
lA = 0.95;
lA_init = lA*ones(1,endNum);
lS_init = lS*ones(1,endNum);
maxiter = 100;
maxiter_ADMM = 50;
method_A = 'ADMM';
psi_est = 1.01;
pos = 1;
norm_A = 1;
param_tau = 0;
display = 0;
S0=param.lib;
endLibScale=cell(1,scaleNum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:scaleNum
    image1=imresize(image1,[size(image2,1),size(image2,2)]);
    x_n1=reshape(image1(:,:,:,i),size(image1,1)*size(image1,2),size(image1,3))';
    x_n2=reshape(image2(:,:,:,i),size(image2,1)*size(image2,2),size(image2,3))';
    %[S,A,TA,psi,scale_factors,J,Rs,Ra,err,t,lA,lS]
%     tmp=hyperFcls(x_n',endLib');
    [S_joint,A_joint,~,psi,scale_factors,~,~,~,~,~,lA,lS] = ...
    bsst({x_n1,x_n2}',S0',size(image2,1),size(image2,2),...
    'maxiter',maxiter,...
    'maxiter_ADMM',maxiter_ADMM,...
    'method_A',method_A,...
    'psi_est',psi_est,...
    'pos',pos,...
    'norm_A',norm_A,...
    'lA_init',lA_init,...
    'lS_init',lS_init,...
    'param_tau',param_tau,...
    'display',display);
    X_hat_tv_i(:,:,:,i)=reshape(A_joint{2}',size(image2,1),size(image2,2),endNum);
    endLibScale{i}=S_joint{2}';
end
param.endLibScale=endLibScale;
end