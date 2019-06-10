%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PR_dec computes the deconvolved and shading-corrected photon-reassigned images        %                                                               %
%                                                                                       %
% Inputs:                                                                               %
%         Ip_est              : estimated patterns                                      %
%         I_PR                : non-deconvolved photon-reassigned image                 %
%         T_incoherent        : incoherent transfer function                            %
%         reg_xi              : regularizer for photon-reassigned image deconvolution   %
%         reg_shading         : regularizer for shading correction                      %
%                                                                                       %
% Outputs:                                                                              %
%         data_projection_dec : deconvolved photon-reassigned image                     %
%         I_obj_corrected     : shading-corrected photon-reassigned image               %
%         H_eff               : effective transfer function                             %
%                                                                                       %
%           Copyright (C) 2016 Li-Hao Yeh                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [data_projection_dec,H_eff,I_obj_corrected] = PR_dec(Ip_est,I_PR,T_incoherent,reg_xi,reg_shading)

F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));

[N,M,Nimg] = size(Ip_est);


% Deconvolution

H_eff = padarray(interp2(T_incoherent.^3),[1,1],'post');
H_eff = H_eff(N/2+1:N*3/2,M/2+1:M*3/2);
H_eff = abs(H_eff/max(abs(H_eff(:))));

reg = reg_xi;

data_projection_dec = real(iF(F(I_PR).*conj(H_eff)./(abs(H_eff).^2+reg)));
data_projection_dec(data_projection_dec<0) = 0;

% Shading correction
Ip_process = zeros(N,M,Nimg);
I_p_mean = mean(Ip_est,3);


for i = 1:Nimg
    Ip_process(:,:,i) = (Ip_est(:,:,i)-I_p_mean).^2;
end

shading_factor = mean(Ip_process,3);

shading_reg = reg_shading*max(shading_factor(:)).^2;

I_obj_corrected = data_projection_dec.*shading_factor./(shading_factor.^2+shading_reg);

end