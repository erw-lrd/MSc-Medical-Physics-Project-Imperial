% SR_COMBINE does a SR reconstruction and
% performs a multi-freq  reconstruction
% This code is based on this paper https://doi.org/10.1002/mrm.28593
% https://github.com/petelally/SuperOffRes/blob/master/src/SR_recon_demo.m

clear; 
close all; 
clc;

%% Load data & dependencies
addpath('external');  % Add path to ifft2c and fft2c functions
load('image_282_data.mat');   % Load demo dataset

% Get dimensions from the loaded data
[N_lin, N_col, N_cha, N_rep] = size(complexFID_all_coils);


%% Super-Resolution Reconstruction

% Pre-allocate space for the final multi-coil image
im_super_res_all_coils = zeros(N_lin * 5, N_col, N_cha, 'like', 1i);

% Loop over all coils
for c = 1:N_cha
    
    phase_step_deg = 72;
    current_coil_kspace = squeeze(complexFID_all_coils(:,:,c,:));
    
    % Initial F-State Separation
    ks_rep0 = current_coil_kspace(:,:,1);
    ks_rep1 = current_coil_kspace(:,:,2);
    ks_rep2 = current_coil_kspace(:,:,3);
    ks_rep3 = current_coil_kspace(:,:,4);
    ks_rep4 = current_coil_kspace(:,:,5);
    
   
    % Calculate F0 State
    ks_F0_initial = (ks_rep0 + ks_rep1 + ks_rep2 + ks_rep3 + ks_rep4) / 5;
    
    % Calculate F-1 State 
    ks_rep1_corr_m1 = ks_rep1 .* exp( 1i * deg2rad(phase_step_deg * 1));
    ks_rep2_corr_m1 = ks_rep2 .* exp( 1i * deg2rad(phase_step_deg * 2));
    ks_rep3_corr_m1 = ks_rep3 .* exp( 1i * deg2rad(phase_step_deg * 3));
    ks_rep4_corr_m1 = ks_rep4 .* exp( 1i * deg2rad(phase_step_deg * 4));
    ks_Fm1_initial = (ks_rep0 + ks_rep1_corr_m1 + ks_rep2_corr_m1 + ks_rep3_corr_m1 + ks_rep4_corr_m1) / 5;
    
    % Calculate F+1 State 
    ks_rep1_corr_p1 = ks_rep1 .* exp(-1i * deg2rad(phase_step_deg * 1));
    ks_rep2_corr_p1 = ks_rep2 .* exp(-1i * deg2rad(phase_step_deg * 2));
    ks_rep3_corr_p1 = ks_rep3 .* exp(-1i * deg2rad(phase_step_deg * 3));
    ks_rep4_corr_p1 = ks_rep4 .* exp(-1i * deg2rad(phase_step_deg * 4));
    ks_F1_initial = (ks_rep0 + ks_rep1_corr_p1 + ks_rep2_corr_p1 + ks_rep3_corr_p1 + ks_rep4_corr_p1) / 5;
    
    % Calculate F-2 State
    ks_rep1_corr_m2 = ks_rep1 .* exp( 2i * deg2rad(phase_step_deg * 1));
    ks_rep2_corr_m2 = ks_rep2 .* exp( 2i * deg2rad(phase_step_deg * 2));
    ks_rep3_corr_m2 = ks_rep3 .* exp( 2i * deg2rad(phase_step_deg * 3));
    ks_rep4_corr_m2 = ks_rep4 .* exp( 2i * deg2rad(phase_step_deg * 4));
    ks_Fm2_initial = (ks_rep0 + ks_rep1_corr_m2 + ks_rep2_corr_m2 + ks_rep3_corr_m2 + ks_rep4_corr_m2) / 5;
    
    % Calculate F+2 State
    ks_rep1_corr_p2 = ks_rep1 .* exp(-2i * deg2rad(phase_step_deg * 1));
    ks_rep2_corr_p2 = ks_rep2 .* exp(-2i * deg2rad(phase_step_deg * 2));
    ks_rep3_corr_p2 = ks_rep3 .* exp(-2i * deg2rad(phase_step_deg * 3));
    ks_rep4_corr_p2 = ks_rep4 .* exp(-2i * deg2rad(phase_step_deg * 4));
    ks_F2_initial = (ks_rep0 + ks_rep1_corr_p2 + ks_rep2_corr_p2 + ks_rep3_corr_p2 + ks_rep4_corr_p2) / 5;


    im_F0_initial  = ifft2c(ks_F0_initial);
    im_F1_initial  = ifft2c(ks_F1_initial);
    im_Fm1_initial = ifft2c(ks_Fm1_initial);
    im_F2_initial  = ifft2c(ks_F2_initial);
    im_Fm2_initial = ifft2c(ks_Fm2_initial);
    
    % Multi-frequency B0 Correction
    % Define phase correction parameters
    num_phase_steps = 100; % Number of phase values to test
    phase_steps = linspace(0, 2*pi, num_phase_steps);

    sr_image_stack = zeros(N_lin * 5, N_col, num_phase_steps, 'like', 1i);
    
    % Loop through each phase step, perform a full SR reco, and store it
    for p = 1:num_phase_steps

        correction_ramp = exp(-1i * phase_steps(p));

        im_F0_corrected  = im_F0_initial  .* correction_ramp;
        im_F1_corrected  = im_F1_initial  .* correction_ramp;
        im_Fm1_corrected = im_Fm1_initial .* correction_ramp;
        im_F2_corrected  = im_F2_initial  .* correction_ramp;
        im_Fm2_corrected = im_Fm2_initial .* correction_ramp;

        ks_F0  = fft2c(im_F0_corrected);
        ks_F1  = fft2c(im_F1_corrected);
        ks_Fm1 = fft2c(im_Fm1_corrected);
        ks_F2  = fft2c(im_F2_corrected);
        ks_Fm2 = fft2c(im_Fm2_corrected);

        % Assembling the SR k-space
        ks_super_res = zeros(N_lin * 5, N_col, 'like', 1i);
        loc_F0  = (N_lin * 2) + 1 : (N_lin * 3);
        loc_F1  = (N_lin * 1) + 1 : (N_lin * 2);
        loc_Fm1 = (N_lin * 3) + 1 : (N_lin * 4);
        loc_F2  = (N_lin * 0) + 1 : (N_lin * 1);
        loc_Fm2 = (N_lin * 4) + 1 : (N_lin * 5);
        
        ks_super_res(loc_F0,  :) = ks_F0;
        ks_super_res(loc_F1,  :) = ks_F1;
        ks_super_res(loc_Fm1, :) = ks_Fm1;
        ks_super_res(loc_F2,  :) = ks_F2;
        ks_super_res(loc_Fm2, :) = ks_Fm2;

        sr_image_stack(:,:,p) = ifft2c(ks_super_res);
    end
    
    % Find the optimal phase step for each pixel in the SR image
    cost_metric = abs(angle(sr_image_stack));
    [~, opt_idx] = min(cost_metric, [], 3);
    
    % Build the final SR image for this coil by picking the best pixels
    im_super_res_single_coil = zeros(N_lin * 5, N_col, 'like', 1i);
    for r = 1:(N_lin * 5)
        for col = 1:N_col
            im_super_res_single_coil(r, col) = sr_image_stack(r, col, opt_idx(r, col));
        end
    end
    
    im_super_res_all_coils(:,:,c) = im_super_res_single_coil;
    
end

% Combine all coil images using Sum of Squares (SoS)
final_im_SoS = sqrt(sum(abs(im_super_res_all_coils).^2, 3));

%% Mean of Low-Resolution Images

mean_k_space = mean(complexFID_all_coils, 4); 
mean_coil_images = zeros(N_lin, N_col, N_cha, 'like', 1i);

for c = 1:N_cha
    k_space_channel = mean_k_space(:, :, c);
    mean_coil_images(:, :, c) = ifft2c(k_space_channel);
end

% Combine coils using SoS
final_image_mean_low_res = sqrt(sum(abs(mean_coil_images).^2, 3));

%% Zero-Filled Image from Repetition 1

% Define the new size for zero-filling (same as SR image)
new_N_lin = N_lin * 5;

k_space_rep1 = complexFID_all_coils(:, :, :, 1);
zf_k_space_all_channels = zeros(new_N_lin, N_col, N_cha, 'like', 1i);

% Place original data
start_row = floor((new_N_lin - N_lin) / 2) + 1;
end_row = start_row + N_lin - 1;
zf_k_space_all_channels(start_row:end_row, :, :) = k_space_rep1;


zf_coil_images = zeros(new_N_lin, N_col, N_cha, 'like', 1i);

for c = 1:N_cha
    k_space_channel_zf = zf_k_space_all_channels(:, :, c);
    zf_coil_images(:, :, c) = ifft2c(k_space_channel_zf);
end

% Combine coils using SoS
final_image_zf = sqrt(sum(abs(zf_coil_images).^2, 3));

%% Display the results

figure('Name', 'Combined Reconstruction Results', 'Position', [100, 100, 1500, 500]);

% Super-Resolution Image
subplot(1, 3, 1);
imagesc(flipud(final_im_SoS));
axis image; 
colormap gray; 
colorbar;
title('Final Super-Resolution Image');

% Mean of Low-Resolution Images
subplot(1, 3, 2);
imagesc(flipud(final_image_mean_low_res));
axis square; 
colormap gray; 
colorbar;
title('Mean of Low-Resolution Images');

% Zero-Filled Image
subplot(1, 3, 3);
imagesc(flipud(final_image_zf));
axis image; 
colormap gray; 
colorbar;
title('Zero-Filled Image (from Rep 1)');

% B0 Phase Offset Plot for the LAST coil processed
figure('Name', 'B0 Phase Offset');
imagesc(flipud(opt_idx));
axis image; 
colormap jet;
colorbar;
title('Optimal B0 phase offset per pixel (Last Coil)');