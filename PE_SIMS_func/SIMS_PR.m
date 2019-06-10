%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMS_PR computes the non-deconvolved photon-reassigned image         %                                                               %
%                                                                      %
% Inputs:                                                              %
%         Ip_est          : estimated patterns                         %
%         I_image_up      : upsampled measurements                     %
%         range_ps        : number of pixel for photon-reassignment    %
%                                                                      %
% Outputs:                                                             %
%         data_projection : shifted covariance image                   %
%         I_PR            : non-deconvolved photon-reassigned image    %
%                                                                      %
%           Copyright (C) 2016 Li-Hao Yeh                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data_projection,I_PR] = SIMS_PR(Ip_est,I_image_up,range_ps)


F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));

[N,M,Nimg] = size(I_image_up);


% Pixel reassignment

% Shift
Nsx = range_ps;
Nsy = range_ps;

Ns = Nsx*Nsy;
pixel_step = 2;

sx = (-(Nsx-1)/2:(Nsx-1)/2).*pixel_step;
sy = (-(Nsy-1)/2:(Nsy-1)/2).*pixel_step;

[syy,sxx] = meshgrid(sy,sx);

s_stack = zeros(2,Ns);

s_stack(1,:) = syy(:);
s_stack(2,:) = sxx(:);



data_projection = zeros(N,M,Nsy,Nsx);

I_p_mean = mean(Ip_est,3);
I_mean = mean(I_image_up,3);

for i = 1:Ns
    [x,y] = ind2sub([Nsx,Nsy],i);
    for j = 1:Nimg
        data_projection(:,:,y,x) = data_projection(:,:,y,x) + circshift((Ip_est(:,:,j)-I_p_mean),[s_stack(1,i),s_stack(2,i)]).*(I_image_up(:,:,j)-I_mean);
    end
end

I_PR = zeros(N,M);

for i = 1:Ns
    [x,y] = ind2sub([Nsx,Nsy],i);
    I_PR = I_PR + circshift(data_projection(:,:,y,x),[-s_stack(1,i)/2,-s_stack(2,i)/2]);
end

end
