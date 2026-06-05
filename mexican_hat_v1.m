%% mexican_hat_v1.m
%  Center-surround ("Mexican hat") receptive-field filtering as a spatial-
%  frequency / edge analyzer -- the visual-system analog of an MRI's k-space.
%
%  A retinal/LGN center-surround cell = excitatory center + inhibitory
%  surround = a Difference-of-Gaussians (DoG), which is a BAND-PASS spatial-
%  frequency filter. Convolving an image with it suppresses flat regions
%  (low spatial frequency) and the finest noise (very high frequency) and
%  keeps edges/texture at the scale set by the center size. Sweeping the
%  center size = tuning to a different spatial-frequency band = the multi-
%  scale channels of early vision.
%
%  Uses base MATLAB only (no Image Processing Toolbox). conv2 does the work.

clear; clc; close all;

%% ---- Settings ----------------------------------------------------------
fname  = 'mriSquare.png';
sigmaC = 2.0;          % center (excitatory) Gaussian width, in pixels
ratio  = 1.6;          % surround/center width ratio (~1.6 approximates a LoG)
scales = [1 2 4 8];    % center widths for the multi-scale channel bank

%% ---- Load image (grayscale double) ------------------------------------
I = double(imread(fname));
if ndims(I) == 3
    I = 0.2989*I(:,:,1) + 0.5870*I(:,:,2) + 0.1140*I(:,:,3);
end
I = I / max(I(:));

%% ---- Build a Mexican-hat (DoG) kernel ----------------------------------
% Helper: a unit-area 2D Gaussian on a grid sized to the surround.
makeDoG = @(sc) buildDoG(sc, ratio);

[DoG, prof, xprof] = makeDoG(sigmaC);

%% ---- Apply it: center-surround / edge response -------------------------
R = conv2(I, DoG, 'same');     % the "V1-ish" response (signed)

%% ---- Figure 1: the kernel and its 1-D profile -------------------------
figure('Color','w','Name','Mexican-hat receptive field','Position',[80 80 900 360]);

subplot(1,3,1);
imagesc(DoG); axis image off; colormap(gca, gray);
title(sprintf('DoG kernel (\\sigma_c=%.1f)', sigmaC));

subplot(1,3,2);
plot(xprof, prof, 'b', 'LineWidth', 1.5); hold on;
yline(0, 'k:'); grid on; axis tight;
xlabel('pixels'); ylabel('weight');
title('1-D profile: + center, - surround');

subplot(1,3,3);
imagesc(R); axis image off; colormap(gca, gray);
title('Center-surround response');

%% ---- Figure 1b: the 3D Mexican hat (sombrero surface) -----------------
% The canonical 3D "sombrero" is the Ricker wavelet (a.k.a. Laplacian-of-
% Gaussian, Marr's grad^2 G) -- the idealized limit of a center-surround
% receptive field, with a tall excitatory center and a pronounced
% inhibitory ring that dips below zero before returning to flat.
sd = sigmaC;  ss = ratio * sigmaC;     % (ss only used in the title)
L  = 4 * sd;  n = 240;
[xg, yg] = meshgrid(linspace(-L, L, n));
r2 = xg.^2 + yg.^2;
H  = (1/(pi*sd^4)) .* (1 - r2/(2*sd^2)) .* exp(-r2/(2*sd^2));   % Ricker / LoG

figure('Color','w','Name','3D Mexican hat (V1 center-surround RF)','Position',[120 120 640 520]);
surf(xg, yg, H, 'EdgeColor','none');
shading interp;
% Diverging blue-white-red colormap centered at zero so SIGN = COLOR:
%   red = excitatory center (+),  blue = inhibitory surround (-).
hk = 128;
blueWhite = [linspace(0.23,1,hk)', linspace(0.33,1,hk)', linspace(0.70,1,hk)'];
whiteRed  = [linspace(1,0.78,hk)', linspace(1,0.10,hk)', linspace(1,0.14,hk)'];
colormap([blueWhite; whiteRed]);
m = max(abs(H(:)));  caxis([-m m]);    % symmetric so white sits exactly at 0
hold on;
% zero plane to make the negative surround ring obvious
zmin = min(H(:));  zmax = max(H(:));
surf(xg, yg, zeros(size(H)), 'FaceAlpha', 0.06, 'EdgeColor','none', 'FaceColor',[.4 .4 .4]);
% zero-crossing contour (boundary between center and surround)
contour3(xg, yg, H, [0 0], 'k-', 'LineWidth', 0.75);
axis tight; grid on;
view(-37, 30);
camlight headlight; lighting gouraud; material dull;
xlabel('x (pixels)'); ylabel('y (pixels)'); zlabel('response  (+ excite / - inhibit)');
title(sprintf('3D Mexican hat: Ricker / Laplacian-of-Gaussian (\\sigma=%.1f)', sd));
% emphasize the negative brim a bit in the vertical view
zlim([1.8*zmin, 1.1*zmax]);

% Save a copy for the book if wanted.
exportName = 'mriSquare_mexhat3D.png';
try
    exportgraphics(gcf, exportName, 'Resolution', 200);   % R2020a+
    fprintf('Saved: %s\n', exportName);
catch
    print(gcf, exportName, '-dpng', '-r200');             % older MATLAB fallback
    fprintf('Saved: %s\n', exportName);
end

%% ---- Figure 2: multi-scale spatial-frequency channels -----------------
figure('Color','w','Name','Multi-scale channels','Position',[100 100 1000 300]);
for i = 1:numel(scales)
    Ki = buildDoG(scales(i), ratio);
    Ri = conv2(I, Ki, 'same');
    subplot(1, numel(scales), i);
    imagesc(Ri); axis image off; colormap(gca, gray);
    title(sprintf('\\sigma_c = %d px', scales(i)));
end
sgtitle('Spatial-frequency channels: small \sigma = fine detail, large \sigma = coarse structure');

%% ---- Save the single-scale edge response ------------------------------
Rdisp = R - min(R(:)); Rdisp = Rdisp / max(Rdisp(:));   % normalize to [0,1] for saving
imwrite(uint8(255*Rdisp), 'mriSquare_mexhat.png');
fprintf('Saved: mriSquare_mexhat.png\n');

%% ---- Local function: build a DoG kernel -------------------------------
function [DoG, prof, xprof] = buildDoG(sigmaC, ratio)
    sigmaS = ratio * sigmaC;                 % surround width
    half   = ceil(3 * sigmaS);               % kernel half-width
    [X, Y] = meshgrid(-half:half, -half:half);
    G  = @(s) exp(-(X.^2 + Y.^2)/(2*s^2));
    Gc = G(sigmaC); Gc = Gc / sum(Gc(:));    % unit-area center
    Gs = G(sigmaS); Gs = Gs / sum(Gs(:));    % unit-area surround
    DoG = Gc - Gs;                           % zero-DC band-pass kernel
    mid  = half + 1;                         % central row for the 1-D profile
    prof = DoG(mid, :);
    xprof = -half:half;
end
