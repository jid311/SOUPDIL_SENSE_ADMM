% Copyright © 2022 by Yan Su and Jizhong Duan, Kunming University of
% Science and Technology. All Rights Reserved.
% This is implementation for SOUPDIL-SENSE

warning off; clc; close all; clear;
addpath(genpath(pwd));

ks = 6;
mc = 24; nc = mc;
dname = 'train_brain_AXT1POST_200_6001959_01_aecc8'; % P3_imk_100 % brain_chall_tr4_tran_95 % image3c_c256 % train_brain_AXT1POST_200_6001959_01_aecc8
load(dname);
DATA = double(DATA);
[m,n,ch] = size(DATA);

matnameroi = sprintf('maskroi_%s', dname); 
load(matnameroi);

for R=[5]

if (m==256)&&(n==256)     load('mask_poisson_256_4.mat');
elseif (m==218)&&(n==170) load('mask_poisson_218x170.mat');
elseif (m==320)&&(n==320) load('mask_poisson_320.mat');
end

mask = mask(:,:,R);
mask((m-mc)/2+1:(m-mc)/2+mc,(n-nc)/2+1:(n-nc)/2+nc) = 1;
sr = nnz(mask)/m/n;

%% Prepare DATA
% DATA = DATA/max(max(sos(ifft2c(DATA))));
DATA = NormalizeImage(DATA);
im = ifft2c(DATA); % get the golden truth image
imr = sos(im);

ksize = [ks,ks]; % ESPIRiT kernel-window-size.
eigThresh_k = 0.02; % threshold of eigenvectors in k-space
eigThresh_im = 0.9; % threshold of eigenvectors in image space

[sx,sy,Nc] = size(DATA);
ncalib = getCalibSize(mask);
mask = repmat(mask,[1,1,Nc]);
DATAc = DATA.*mask;
calib = crop(DATAc,[ncalib,Nc]);

%% Compute Eigen-Value Maps
% Maps are computed in two steps. Compute Calibration matrix, perform 1st SVD and convert singular vectors into k-space kernels
[k,S] = dat2Kernel(calib,ksize);
idx = max(find(S >= S(1)*eigThresh_k));

%% crop kernels and compute eigen-value decomposition in image space to get maps
[M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

%% Compute Soft-SENSE ESPIRiT Maps
% crop sensitivity maps according to eigenvalues==1. Note that we have to % use 2 sets of maps. Here we weight the 2 maps with the eigen-values
maps = M(:,:,:,end-1:end);

% Weight the eigenvectors with soft-senses eigen-values
weights = W(:,:,end-1:end) ;
weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end-1:end) > eigThresh_im);
weights = -cos(pi*weights)/2 + 1/2;

% create and ESPIRiT operator
ESP = ESPIRiT(maps,weights);
SNS = ESPIRiT(maps(:,:,:,end));
FT = p2DFT(mask,[sx,sy,Nc]);

resSENSE = zeros(sx,sy);
x0 = SNS'*(FT'.*DATAc);

opts = struct('FT',FT,'DATA',DATA,'maps',maps(:,:,:,end),'mask',mask,'maskroi',maskroi,'im',imr, ...
'mu1',1,'mu2',1,'a2',1,'x0',x0,'num2',1,'wd',1,'n',36,'IterOG',70);

% Parameters should be tuned for different datasets and accelerating facters (R)
tab_a = [0.3]; tab_a2 = [0.058]; tab_mu1 = [0.12]; tab_mu2 = [0.01]; opts.IterOG = 4; % P3_imk_100, R=5
% tab_a = [0.12]; tab_a2 = [0.0422]; tab_mu1 = [0.086]; tab_mu2 = [0.01]; opts.IterOG = 100; % brain_chall_tr4_tran_95, R=5
% tab_a = [0.6]; tab_a2 = [0.024]; tab_mu1 = [0.079]; tab_mu2 = [0.013]; opts.IterOG = 180; % image3c_c256, R=5  
% tab_a = [0.13]; tab_a2 = [0.03]; tab_mu1 = [0.2]; tab_mu2 = [0.013]; opts.IterOG = 100; % train_brain_AXT1POST_200_6001959_01_aecc8, R=5

for imu1 = 1:numel(tab_mu1)
for imu2 = 1:numel(tab_mu2)
for ia = 1:numel(tab_a)
for ia2 = 1:numel(tab_a2)

opts.mu1 = tab_mu1(imu1);
opts.mu2 = tab_mu2(imu2);
opts.a = tab_a(ia);
opts.a2 = tab_a2(ia2);

[res2, out] = SOUPDIL_SENSE_ADMM(DATAc, opts);
plot(out.snr);

fprintf('%s, R=%d, IterOG:%g, a:%g, a2:%g, mu1:%g, mu2:%g, MSNR:%.2f, ESNR:%.2f, EHFEN:%.4f, ESSIM:%.4f, time: %.0f\n',...
dname,R, opts.IterOG, opts.a, opts.a2, opts.mu1,opts.mu2,max(out.snr),out.snr(end),out.hfen(end),out.ssim(end), out.timet(end));
end
end
end
end

end
