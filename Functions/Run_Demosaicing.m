function [I_HS] = Run_Demosaicing(I,num_band,Algo)
%% This function is used to apply demosaicing to an input SSI image
% Author : Kinan ABBAS
% Creation Date : 3 Mars 2023
%collaborators:DAHOU SOUKAYNA,OUMENSKOU YOUSSEF,AZOUAOUI MERIEME

% Input: 
%    I : The 2D SSI image
%    num_bands: number of bands in the image. The possible values are 16
%    and 25
%    Algo: the prefered demosaicing algorithm. The possible values are : 1
%    for Weighted Bilinear Interpolation

% Output:
%   I_HS : The restored 3D datacube


% Scale the values to be between 1 and 255, this step can be ignored
I=I/255;


% Get the size of the image
[n1,n2,n3]=size(I);

% Set the algorithm to default value
if(isempty(Algo))
    Algo=1;
end

% Protecting step to ensure the image size is correct
if num_band==25
    if(mod(n1,5)~=0 || mod(n2,5)~=0)
        disp('Error! SSI image dimensions must be divided by 5 !')
    end
elseif num_band==16
    if(mod(n1,4)~=0 || mod(n2,4)~=0)
        disp('Error! SSI image dimensions must be divided by 5 !')
    end
else
    disp('Error ! choose the correct number of bands');
end

% Creat the filter pattern
[~,FilterPattern_lst]=make_sampling_operators2(n1,n2,num_band,1,num_band,1);
FilterPattern=cell2mat(FilterPattern_lst);

% Run demosaicing
switch Algo
    case 1
        
        I_WB_tmp=run_WB(I,FilterPattern,num_band);
        I_HS=mean(I_WB_tmp,4);
    case 2
        disp('Running WNMF')
        % To add the implementation of WNMF here
    case 3
        disp('Running PPID');
        %I_PPID_tmp=zeros(sz(1),sz(2),num_band,num_obs_pxl);
        % To add the implementation of PPID here
        
        
    otherwise
        disp('The selected algorithm can not be found');
end
    


end

