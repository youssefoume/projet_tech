% %% A Script to apply demosaicing to an SSI image
% % Author: Kinan ABBAS
% % Creation Date: 3 Mars 2023
% %collaborators:DAHOU SOUKAYNA,OUMENSKOU YOUSSEF,AZOUAOUI MERIEME
% 
% 
% 
% % Cleaning and loading everthing to the path
% close all;
% clear;
% 
% addpath(pwd);
% cd Data/
% addpath(genpath(pwd));
% cd ..;

cd Functions/
addpath(genpath(pwd));
cd ..;

cd Methods/
addpath(genpath(pwd));
cd ..;

% 
% %% Load the SSI image
image = imread("U:\BD_TECH\cup\HSI\0692.png"); % The function will load the SSI image data. The image is stored in the 'image' variable
% 
% % Show the image
% figure; imagesc(image);title('The SSI image'); colorbar;
% 
% % Convert image to double
% I=double(image);
% 
% 
% %% Demosaicing using WB
% % Call the demosaicing function with WB method
% I_HS_WB= Run_Demosaicing(I(:,:,1),16,1);
% 
% % Show the first two bands of the restored datacube
% figure;
% subplot(1,2,1); imagesc(I_HS_WB(:,:,1));title('Band 1');
% subplot(1,2,2); imagesc(I_HS_WB(:,:,2)); title('Band 2');
% 
% 
I_HS_WB=Demosaicing_img(image);
figure;
subplot(1,2,1); imagesc(I_HS_WB(:,:,1));title('Band 1');
subplot(1,2,2); imagesc(I_HS_WB(:,:,2)); title('Band 2');
