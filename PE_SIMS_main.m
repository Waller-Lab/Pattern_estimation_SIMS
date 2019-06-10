%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Pattern estimation SIM with statistical prior    %
%           Copyright (C) 2017 Li-Hao Yeh           %
%                                                   %
% Please cite:                                      %
%                                                   %
% L.-H. Yeh, L. Tian, and L. Waller, "Structured    %
% illumination microscopy with unknown patterns and %
% a statistical prior," Biomed. Opt. Express 8,     %
% 695-711 (2017)                                    %
%                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
set(0,'DefaultFigureWindowStyle','docked');

addpath('PE_SIMS_func');



F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));

load res_speckle_shift.mat;  % Simulated data

% Microtubule data from A. York et al., Nat. Methods 2013.
% load Nat_tubule_data.mat;  

%% Coordinate assignment

[Ncrop,Mcrop,Nimg] = size(I_image);


N = Ncrop*2; M = Mcrop*2; ps = pscrop/2;

xh = (-M/2:(M/2-1)).*ps; yh = (-N/2:(N/2-1)).*ps;
fx = (-M/2:(M/2-1))./(ps*M); fy = (-N/2:(N/2-1))./(ps*N);
NAx = fx*lambda; NAy = fy*lambda;
[xhh,yhh] = meshgrid(xh,yh);
[fxx,fyy] = meshgrid(fx,fy);

%% Upsampling the data

I_image_up = zeros(N,M,Nimg);
for i = 1:Nimg
    temp = max(0,I_image(:,:,i) - 100);
    I_image_up(:,:,i) = abs(iF(padarray(F(temp),[(N-Ncrop)/2,(M-Mcrop)/2])));
end

%% Iterative Pattern Estimation

% incoherent transfer function

Pupil_obj = zeros(N,M);
r_obj=(fxx.^2+fyy.^2).^(1/2);
Pupil_obj(find(r_obj<NA_obj/lambda))=1;
T_incoherent = abs(F(abs(iF(Pupil_obj)).^2));
T_incoherent = T_incoherent/max(T_incoherent(:));

%% Pattern estimation

reg_beta = 1e-3;
reg_delta = 1e-3;
max_itr = 50;

[Ip_est,I_mean,I_mean_dec] = IPE(I_image_up,T_incoherent,max_itr,reg_beta,reg_delta);

%% SIMS process


reg_xi = 1e-4;
reg_shading = 1e-2;
[data_projection,data_projection_dec,I_obj_corrected,H_eff] = SIMS(Ip_est,I_image_up,T_incoherent,reg_xi,reg_shading);

%% SIMS-PR process

if mod(ceil(lambda/2/NA_obj/ps),2) == 0 
    range_ps = ceil(lambda/2/NA_obj/ps)+1;
else
    range_ps = ceil(lambda/2/NA_obj/ps);
end
[data_projection_pr,I_PR] = SIMS_PR(Ip_est,I_image_up,range_ps);

%% SIMS-PR deconvolution and shading correction

reg_xi_pr = 1e-5;
reg_shading_pr = 1e-2;

[I_PR_dec,H_eff_PR,I_PR_corrected] = PR_dec(Ip_est,I_PR,T_incoherent,reg_xi_pr,reg_shading_pr);

%% Plot result

figure;imagesc(I_mean_dec);colormap gray;axis square;
title('deconvolved mean raw data');
figure;imagesc(data_projection_dec);colormap gray;axis image;
title('IPE SIMS reconstruction');
figure;imagesc(I_obj_corrected);colormap gray;axis image;
title('IPE SIMS reconstruction with shading correction');

figure;imagesc(I_PR_dec);colormap gray;axis image;
title('IPE SIMS-PR reconstruction');
figure;imagesc(I_PR_corrected);colormap gray;axis image;
title('IPE SIMS-PR reconstruction with shading correction');
