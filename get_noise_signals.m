function xn = get_noise_signals (COV, Nsamples)

% GET_NOISE_SIGNALS: Generates noise signals from a noise covariance matrix
% INPUT:
% - COV: Noise covariance matrix (M x M)
% - Nsamples: Number of time points (length of noise signals)
% OUTPUT:
% - xn: noise signals (M x Nsamples)
%
% Author: Guiomar Niso, 2014
%

[V,D] = eig(COV);

% xn = (1/SNR) * V * D.^(1/2) * randn(size(COV,1),Nsamples);
xn = V * D.^(1/2) * randn(size(COV,1),Nsamples);


%%%%%%
% Example:
% SNR = 0.3;
% Nsamples = 500;
% xn = get_noise_signals (n.NoiseCov, Nsamples);
% xnn = xn./max(max(xn));
% xns = xnn.*max(max(s.F));
% sn = s.F + SNR*xns;
% s.F=sn;

% figure(1); imagesc(n.NoiseCov); colorbar;
% figure(2); imagesc(xn*xn' ./ (size(xn,2) - 1)); colorbar;
% figure(3); imagesc(cov(xn)); colorbar;
% See also noise extracted from recordings
