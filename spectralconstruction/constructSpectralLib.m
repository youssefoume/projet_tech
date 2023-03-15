function [lib,endmembers]=constructSpectralLib(img_path,dicNum,block,bands)
% construct spectral lib from video
% img_path: img files in png format
% dicNum: the spectral num
img_files = dir(fullfile(img_path, '*.png'));
img_files = {img_files.name};
img_files=cellstr(img_files);
imgNum=length(img_files);
verbose='off';
endmembers=cell(1,imgNum);
for i=1:imgNum
    im = imread([img_path '/' img_files{i}]);
    im=X2Cube(double(im),block,bands);%convert PNG to hsi
    im=reshape(im,size(im,1)*size(im,2),size(im,3));
    kf=estimateR(im','additive',verbose);
    endmembers{i}=hyperVca(im',kf);
end
endmembers=cell2mat(endmembers);
[idex,lib]=kmeans(endmembers, dicNum);
lib=lib';
end