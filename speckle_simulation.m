%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Simulated data generation for PE-SIMS       %
%           Copyright (C) 2016 Li-Hao Yeh           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
set(0,'DefaultFigureWindowStyle','docked');



F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));

I = double(imread('resolution.jpg'));


I = I(323:506,358:541,2);
I = I/max(I(:));



%% --------Experiment Setup----------

lambda = 0.605; k=2*pi/lambda; % wavelength (micons) and wave number
mag = 120;
pscrop = 6.5/mag; % Pixels size (microns)
NA_obj = 0.8;
NAs = 0.8;

undersamp_factor = 1;
ps = pscrop/undersamp_factor;


[N,M] = size(I);

Ncrop = N/undersamp_factor;
Mcrop = M/undersamp_factor;

xh = (-M/2:(M/2-1)).*ps; yh = (-N/2:(N/2-1)).*ps;
fx = (-M/2:(M/2-1))./(ps*M); fy = (-N/2:(N/2-1))./(ps*N);
NAx = fx*lambda; NAy = fy*lambda;
[xhh,yhh] = meshgrid(xh,yh);
[fxx,fyy] = meshgrid(fx,fy);


%% Star generator
% [theta,rho] = cart2pol(xhh(:),yhh(:));
% 
% theta = reshape(theta,N,M);
% rho = reshape(rho,N,M);
% 
% I = 1+cos(40*theta);
% I = I(11:N-10,11:M-10);
% I = padarray(I,[10,10]);
% 
% figure;imagesc(I);colormap gray;axis image;axis off;

%% -------- Propagation kernel --------

xs = (-M*4/2:(M*4/2-1)).*ps; ys = (-N*4/2:(N*4/2-1)).*ps;
fxs = (-M*4/2:(M*4/2-1))./(ps*4*M); fys = (-N*4/2:(N*4/2-1))./(ps*4*N);
NAxs = fxs.*lambda; NAys = fys.*lambda;
[xss,yss] = meshgrid(xs,ys);
[fxxs,fyys] = meshgrid(fxs,fys);
NAxx = fxxs.*lambda;
NAyy = fyys.*lambda;


r_prop=(fxxs.^2+fyys.^2).^(1/2);
Pupil_prop = zeros(N*4,M*4);
Pupil_prop(find(r_prop<NAs/lambda))=1;


Pupil_2NA = zeros(N*4,M*4);
Pupil_2NA(find(r_prop<2*NAs/lambda))=1;


%% Pattern generation with random phase mask


rng(50);

pattern_idx = 1; % 1: multi-spot, 2: speckle_pattern


switch pattern_idx
    case 1
        dot_period = 24;
        window_size = 2;
        % 
        speckle_intensity = 0.01*ones(4*N,4*M);

        num_period_y = floor(4*N/dot_period);
        num_period_x = floor(4*M/dot_period);

        for i = 1:num_period_y
            for j = 1:num_period_x
                idx_y = i*dot_period-dot_period/2;
                idx_x = j*dot_period - dot_period/2;
                speckle_intensity((idx_y-(window_size/2-1)):(idx_y+(window_size/2)),(idx_x-(window_size/2-1)):(idx_x+(window_size/2))) = 1;
            end
        end

        % 
        speckle_intensityf = F(speckle_intensity);
        speckle_intensity = abs(iF(speckle_intensityf.*Pupil_prop)).^2;
        
    case 2
        random_map = rand(4*N,4*M);

        random_map(random_map<=0.9) = 0;
        random_map(random_map>0.9) = 1;
        % 
        random_mapf = F(random_map);
%         random_mapf(2*N+1,2*M+1)=0; % uncomment if it's fully-developed
        random_mapf = exp(1j*rand(4*N,4*M)*100);
        speckle_intensity = abs(iF(random_mapf.*Pupil_prop)).^2;
end

speckle_intensity = abs(iF(F(speckle_intensity).*Pupil_2NA));
speckle_intensity = speckle_intensity/max(speckle_intensity(:));
speckle_intensity_crop = speckle_intensity(1.5*N+1:2.5*N,1.5*M+1:2.5*M);
speckle_intensity_cropf = F(speckle_intensity_crop);


figure;imagesc(xh,yh,speckle_intensity_crop); colormap gray;axis square;
% % figure;imagesc(NAx,NAy,log10(abs(F(speckle_field_crop))),[1 5]);colormap jet;axis square;
figure;imagesc(NAx,NAy,log10(abs(speckle_intensity_cropf))/max(max(log10(abs(speckle_intensity_cropf)))),[0 1]); colormap jet; axis square;
hold on;circle(0,0,2*NAs);


%% Generate simulated data with angular/pixel shift

% Pixel shift part

N_shiftx = 6;
N_shifty = 6;

Nimg = N_shiftx*N_shifty;
pixel_step = 4;

pixel_shiftx = (-(N_shiftx-1)/2:(N_shiftx-1)/2).*pixel_step;
pixel_shifty = (-(N_shifty-1)/2:(N_shifty-1)/2).*pixel_step;

[pixel_shiftyy,pixel_shiftxx] = meshgrid(pixel_shifty,pixel_shiftx);

pixel_shift_stack = zeros(2,Nimg);

pixel_shift_stack(1,:) = pixel_shiftyy(:);
pixel_shift_stack(2,:) = pixel_shiftxx(:);

I_image = zeros(Ncrop,Mcrop,Nimg);


%%

Pupil_obj = zeros(N,M);
r_obj=(fxx.^2+fyy.^2).^(1/2);
Pupil_obj(find(r_obj<NA_obj/lambda))=1;
T_incoherent = abs(F(abs(iF(Pupil_obj)).^2));
T_incoherent = T_incoherent/max(T_incoherent(:));

speckle_intensity_shift_crop = zeros(M,N,Nimg);
speckle_intensity_r_shift_crop = zeros(M,N,Nimg);

for i = 1:Nimg
    
    % Pixel shift part
    
    speckle_intensity_shift_crop(:,:,i) = speckle_intensity(1.5*N+1+pixel_shift_stack(1,i):2.5*N+pixel_shift_stack(1,i),1.5*M+1+pixel_shift_stack(2,i):2.5*M+pixel_shift_stack(2,i));
   
    Itemp = F(speckle_intensity_shift_crop(:,:,i).*I).*T_incoherent;
    I_image(:,:,i) = abs(iF(Itemp(N/2+1-Ncrop/2:N/2+1+Ncrop/2-1,M/2+1-Mcrop/2:M/2+1+Mcrop/2-1)));

end

%% Poisson + background noise

dark_current = 100;
photon_count = 5000;

% I_image = I_image/max(I_image(:))*photon_count;

I_image = imnoise(( I_image/max(I_image(:)).*photon_count + dark_current)*1e-12,'poisson')*1e12;


%% Save file
savefile='res_speckle_shift';
save(savefile,'pscrop','lambda','NA_obj','I_image','I','speckle_intensity','pixel_shift_stack','speckle_intensity_shift_crop');
