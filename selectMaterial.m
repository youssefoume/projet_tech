function endLib=selectMaterial(image,endLib)
x_n=reshape(image(:,:,:),size(image,1)*size(image,2),size(image,3));
% lambda = 5e-3;
% lambda_TV = 3e-3;
% endLib=max(endLib,eps);
lambda = 1e-3;
R=estimateR(x_n','additive','off');
% [tmp] =  sunsal(endLib,x_n','lambda',lambda,'ADDONE','yes','POSITIVITY','yes', ...
%                     'TOL',1e-4, 'AL_iters',2000,'verbose','no');
[tmp] = clsunsal(endLib,x_n','POSITIVITY','yes','VERBOSE','no','ADDONE','yes', ...
    'lambda', 3e-4,'AL_ITERS',200, 'TOL', 1e-6);
tmp(tmp<0)=0;
[val,pos]=sort(sum(tmp,2));
endLib=endLib';
endLib= endLib(pos(end-floor(R)+1:end),:);
% endLib=max(endLib,eps);
% endLib=hyperVca(x_n',R)';
% endLib=max(endLib,eps);
% abundance=hyperFcls(x_n',endLib');
% [endLib,abundance]=nnschalfprior(x_n',endLib',abundance);
% endLib=endLib';

end