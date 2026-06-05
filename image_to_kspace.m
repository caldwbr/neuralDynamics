%% image_to_kspace.m
%  Convert an image to its MRI-style k-space (2D spatial-frequency) representation.
%
%  k-space is just the 2D Fourier transform of the image: each point is a
%  spatial-frequency component (the center = low frequencies / overall
%  brightness and contrast; the edges = high frequencies / fine detail and
%  sharp edges). Raw MRI scanners actually acquire this k-space data and then
%  inverse-Fourier-transform it to form the picture; here we go the other way.
%
%  Uses base MATLAB only (no Image Processing Toolbox required).
%  Tested convention: centered FFT so DC (zero frequency) sits in the middle.
%
%  Author: (your name)

clear; clc; close all;

%% ---- 1. Load the image -------------------------------------------------
fname = 'mriSquare.png';
I = imread(fname);

% Convert to double precision grayscale in [0,1].
I = double(I);
if ndims(I) == 3                       % RGB -> luminance (Rec. 601 weights)
    I = 0.2989*I(:,:,1) + 0.5870*I(:,:,2) + 0.1140*I(:,:,3);
end
I = I / max(I(:));                     % normalize (FFT is scale-agnostic, but tidy)

%% ---- 2. Forward transform: image -> k-space ---------------------------
% ifftshift before fft2 and fftshift after keeps the image centered and makes
% the round trip exact. K is complex: it has both magnitude and phase.
K = fftshift(fft2(ifftshift(I)));

mag   = abs(K);                        % k-space magnitude
phase = angle(K);                      % k-space phase (radians)

% Magnitude spans many orders of magnitude, so display it on a log scale.
logMag = log(1 + mag);
logMagDisp = (logMag - min(logMag(:))) / (max(logMag(:)) - min(logMag(:)));

%% ---- 3. Verify: k-space -> image (inverse transform) ------------------
Irecon = real(fftshift(ifft2(ifftshift(K))));
reconError = max(abs(Irecon(:) - I(:)));
fprintf('Round-trip max reconstruction error: %.3e (should be ~1e-15)\n', reconError);

%% ---- 4. Display --------------------------------------------------------
figure('Color','w','Name','Image <-> k-space','Position',[100 100 1100 380]);

subplot(1,3,1);
imagesc(I); axis image off; colormap(gca, gray);
title('Image (spatial domain)');

subplot(1,3,2);
imagesc(logMagDisp); axis image off; colormap(gca, gray);
title('k-space magnitude (log)');

subplot(1,3,3);
imagesc(phase); axis image off; colormap(gca, gray);
title('k-space phase');

%% ---- 5. Save outputs ---------------------------------------------------
imwrite(uint8(255*logMagDisp), 'mriSquare_kspace_mag.png');   % log-magnitude
phaseDisp = (phase + pi) / (2*pi);                            % map [-pi,pi]->[0,1]
imwrite(uint8(255*phaseDisp), 'mriSquare_kspace_phase.png');  % phase

% Save the raw complex k-space so you can reload and reconstruct later.
save('mriSquare_kspace.mat', 'K');

fprintf('Saved: mriSquare_kspace_mag.png, mriSquare_kspace_phase.png, mriSquare_kspace.mat\n');

% % ======================================================================
%  OPTIONAL MRI DEMOS (uncomment to run) — these show *why* k-space matters
%  ======================================================================

% (a) Keep only the CENTER of k-space (low frequencies) -> blurry image.
%     This is what "low-resolution" / undersampled MRI acquisition looks like.
[ny, nx] = size(K);
r = 0.12;                                   % fraction of k-space radius to keep
[xx, yy] = meshgrid(1:nx, 1:ny);
cx = nx/2 + 0.5; cy = ny/2 + 0.5;
lowpass = sqrt(((xx-cx)/(nx/2)).^2 + ((yy-cy)/(ny/2)).^2) <= r;
Iblur = real(fftshift(ifft2(ifftshift(K .* lowpass))));
figure('Color','w'); imagesc(Iblur); axis image off; colormap gray;
title(sprintf('Center %d%% of k-space only (low-pass)', round(100*r)));

% (b) Remove the center (high frequencies only) -> edges/outlines.
Iedge = real(fftshift(ifft2(ifftshift(K .* ~lowpass))));
figure('Color','w'); imagesc(Iedge); axis image off; colormap gray;
title('High-pass k-space (edges only)');
