%% gabor_wavelets.m
%  Gabor wavelets at two spatial frequencies (color-coded), plus an
%  even/odd "quadrature pair" at one frequency.
%
%  A 2D Gabor = a sinusoidal carrier under a Gaussian envelope. It is the
%  standard model of a V1 simple-cell receptive field, tuned to an
%  orientation, a spatial frequency (1/lambda), and a PHASE (phi: 0 = even
%  / cosine, pi/2 = odd / sine). Here a low-frequency Gabor is tinted green
%  and a high-frequency one red; the overlay shows their interference.
%
%  Uses base MATLAB only (no Image Processing Toolbox).

clear; clc; close all;

%% ---- Grid and parameters ----------------------------------------------
N     = 257;
[x,y] = meshgrid(linspace(-1,1,N));
sigma = 0.30;        % Gaussian envelope width
theta = 30;          % orientation (degrees), same for both so we compare frequency
lam1  = 0.45;        % LOW spatial frequency  (long wavelength) -> green
lam2  = 0.18;        % HIGH spatial frequency (short wavelength) -> red

% 2D Gabor: rotate coords, Gaussian envelope * cosine carrier with phase phi.
gabor = @(th,lam,phi) ...
    exp(-(((x*cosd(th)+y*sind(th)).^2) + ((-x*sind(th)+y*cosd(th)).^2))/(2*sigma^2)) ...
    .* cos(2*pi*(x*cosd(th)+y*sind(th))/lam + phi);

nrm = @(g) (g - min(g(:))) / (max(g(:)) - min(g(:)));   % scale to [0,1] for display

%% ---- Two frequencies, color-coded -------------------------------------
g1 = nrm(gabor(theta, lam1, 0));     % low frequency  (even phase)
g2 = nrm(gabor(theta, lam2, 0));     % high frequency (even phase)
z  = zeros(N);

greenImg = cat(3, z,  g1, z);        % low frequency in the GREEN channel
redImg   = cat(3, g2, z,  z);        % high frequency in the RED channel
overlay  = cat(3, g2, g1, z);        % R = high, G = low; overlap -> yellow

figure('Color','k','Name','Gabor wavelets','Position',[100 100 1200 420]);
subplot(1,3,1); image(greenImg); axis image off;
title(sprintf('low frequency  (\\lambda=%.2f)', lam1), 'Color','w');
subplot(1,3,2); image(redImg);   axis image off;
title(sprintf('high frequency (\\lambda=%.2f)', lam2), 'Color','w');
subplot(1,3,3); image(overlay);  axis image off;
title('overlay (interference)', 'Color','w');

%% ---- Even / odd quadrature pair at one frequency ----------------------
% V1 simple cells come in even (cosine) and odd (sine) pairs 90 deg apart in
% phase; a complex cell squares and sums them for a PHASE-INVARIANT response.
ge = gabor(theta, lam1, 0);          % even  (cosine phase)
go = gabor(theta, lam1, pi/2);       % odd   (sine phase)
energy = sqrt(ge.^2 + go.^2);        % quadrature energy = phase-invariant envelope

figure('Color','w','Name','Quadrature pair','Position',[120 120 1100 320]);
subplot(1,3,1); imagesc(ge); axis image off; colormap(gca,gray);
title('even  (\phi = 0)');
subplot(1,3,2); imagesc(go); axis image off; colormap(gca,gray);
title('odd  (\phi = \pi/2)');
subplot(1,3,3); imagesc(energy); axis image off; colormap(gca,gray);
title('quadrature energy (phase-invariant)');

%% ---- Save ---------------------------------------------------------------
imwrite(uint8(255*greenImg), 'gabor_green.png');
imwrite(uint8(255*redImg),   'gabor_red.png');
imwrite(uint8(255*overlay),  'gabor_overlay.png');
fprintf('Saved: gabor_green.png, gabor_red.png, gabor_overlay.png\n');
