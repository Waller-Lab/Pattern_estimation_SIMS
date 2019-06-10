%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMS computes the deconvolved and shading-corrected covariance images          %                                                               %
%                                                                                %
% Inputs:                                                                        %
%         Ip_est              : estimated patterns                               %
%         I_image_up          : upsampled measurements                           %
%         T_incoherent        : incoherent transfer function                     %
%         reg_xi              : regularizer for covariance image deconvolution   %
%         reg_shading         : regularizer for shading correction               %
%                                                                                %
% Outputs:                                                                       %
%         data_projection     : non-deconvolved covariance image                 %
%         data_projection_dec : deconvolved covariance image                     %
%         I_obj_corrected     : shading-corrected covariance image               %
%         H_eff               : effective transfer function                      %
%                                                                                %
%           Copyright (C) 2016 Li-Hao Yeh                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data_projection,data_projection_dec,I_obj_corrected,H_eff] = SIMS(Ip_est,I_image_up,T_incoherent,reg_xi,reg_shading)

F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));

[N,M,Nimg] = size(I_image_up);


data_projection = zeros(N,M);

I_p_mean = mean(Ip_est,3);
I_mean = mean(I_image_up,3);
for j = 1:Nimg
    data_projection = data_projection + (Ip_est(:,:,j)-I_p_mean).*(I_image_up(:,:,j)-I_mean);
end

H_eff = F(iF(T_incoherent.^2).*iF(T_incoherent));
H_eff = abs(H_eff/max(abs(H_eff(:))));

data_projection_dec = iF(F(data_projection).*conj(H_eff)./(abs(H_eff).^2+reg_xi));
data_projection_dec(data_projection_dec<0) = 0;

% Shading correction
Ip_process = zeros(N,M,Nimg);

for i = 1:Nimg
    Ip_process(:,:,i) = (Ip_est(:,:,i)-I_p_mean).^2;
end

shading_factor = mean(Ip_process,3);

shading_reg = reg_shading*max(shading_factor(:))^2;

I_obj_corrected = data_projection_dec.*shading_factor./(shading_factor.^2+shading_reg);
end