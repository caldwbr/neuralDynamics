%% Jordan (1986) square limit-cycle attractor with inside/outside spirals
%
%  Recreates the dynamical picture from Michael Jordan,
%  "Serial Order: A Parallel Distributed Processing Approach" (1986),
%  Figures 8 and 9: a network that has learned a 4-corner cyclic
%  sequence on the unit square. The four corners
%      C1 = (.25,.25)  C2 = (.75,.25)
%      C3 = (.75,.75)  C4 = (.25,.75)
%  form a limit cycle. Trajectories starting inside the square spiral
%  outward to it; trajectories starting outside spiral inward to it.
%
%  Dynamics (discrete, one update per corner-visit):
%     P_{t+1} = C_{next} + gamma * R_90 * (P_t - C_{curr})
%
%  where R_90 is the 90-degree CCW rotation matrix and gamma in (0,1)
%  is the convergence rate. The offset (P_t - C_curr) gets rotated by
%  90 degrees and scaled by gamma each step. Connecting successive
%  P_t with line segments produces the characteristic rectangular
%  spirals Jordan documented.

clear; clc; close all;

%% Geometry
corners = [0.25 0.25;
           0.75 0.25;
           0.75 0.75;
           0.25 0.75];      % cyclic, CCW

R = [0 -1; 1 0];            % 90-degree CCW rotation

%% Trajectory parameters
gamma  = 0.50;              % convergence rate (0 < gamma < 1)
nLaps  = 6;                 % laps of the cycle to simulate
nSteps = 4 * nLaps;

%% Simulate
P_in  = simulate([0.40, 0.40], corners, R, gamma, nSteps);
P_out = simulate([0.05, 0.05], corners, R, gamma, nSteps);

%% Plot
fig = figure('Color','white', 'Position', [80 60 900 820]);
hold on; axis equal; box on; grid on;
xlim([-0.02, 1.02]); ylim([-0.02, 1.02]);
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1, ...
         'FontSize', 12, 'LineWidth', 1.0, ...
         'GridAlpha', 0.25);

% Square attractor (the learned limit cycle)
sq = [corners; corners(1, :)];
hSq = plot(sq(:,1), sq(:,2), '-', 'Color', [0 0 0], 'LineWidth', 3);

% Inside trajectory (spirals outward to the square)
hIn = plot(P_in(:,1), P_in(:,2), '-', ...
           'Color', [0.20 0.40 0.85], 'LineWidth', 1.6);
plot(P_in(1,1), P_in(1,2), 'o', ...
     'MarkerSize', 11, 'LineWidth', 1.2, ...
     'MarkerFaceColor', [0.20 0.40 0.85], 'MarkerEdgeColor', 'k');

% Outside trajectory (spirals inward to the square)
hOut = plot(P_out(:,1), P_out(:,2), '-', ...
            'Color', [0.85 0.30 0.20], 'LineWidth', 1.6);
plot(P_out(1,1), P_out(1,2), 's', ...
     'MarkerSize', 12, 'LineWidth', 1.2, ...
     'MarkerFaceColor', [0.85 0.30 0.20], 'MarkerEdgeColor', 'k');

% Annotate start points
text(0.40, 0.40, '  (.4, .4)', 'FontSize', 11, ...
     'Color', [0.15 0.30 0.65], 'VerticalAlignment','bottom');
text(0.05, 0.05, '  (.05, .05)', 'FontSize', 11, ...
     'Color', [0.65 0.20 0.15], 'VerticalAlignment','bottom');

% Annotate corner labels
cornerLabels = {'(.25, .25)', '(.75, .25)', '(.75, .75)', '(.25, .75)'};
cornerOffsets = [-0.02 -0.04;  0.02 -0.04;  0.02  0.04; -0.02  0.04];
for k = 1:4
    plot(corners(k,1), corners(k,2), 'ko', ...
         'MarkerSize', 6, 'MarkerFaceColor', 'k');
    text(corners(k,1) + cornerOffsets(k,1), ...
         corners(k,2) + cornerOffsets(k,2), ...
         cornerLabels{k}, 'FontSize', 10, ...
         'HorizontalAlignment', 'center');
end

xlabel('unit 1', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('unit 2', 'FontSize', 14, 'FontWeight', 'bold');
title({'Square attractor (limit cycle) with spiraling trajectories', ...
       'after Jordan (1986), Figures 8--9'}, ...
      'FontSize', 13, 'FontWeight', 'bold');

legend([hSq, hIn, hOut], ...
       {'Limit cycle (learned 4-corner sequence)', ...
        'Inside start: spirals outward', ...
        'Outside start: spirals inward'}, ...
        'Location', 'south', 'FontSize', 11);

%% Save
outDir = fileparts(mfilename('fullpath'));
if isempty(outDir), outDir = pwd; end
outFile = fullfile(outDir, 'jordanAttractor.png');
exportgraphics(fig, outFile, 'Resolution', 300);
fprintf('Saved figure to: %s\n', outFile);

%% ================== TRAJECTORY FUNCTION ==================
function P = simulate(startPt, corners, R, gamma, nSteps)
    % Discrete dynamics with one update per corner-visit:
    %   At step t, the trajectory was "near" corner C_curr; the
    %   network now drives it toward C_next, with the residual
    %   offset (P - C_curr) rotated 90 degrees CCW and scaled by gamma.
    %   Connecting consecutive P_t with line segments produces the
    %   rectangular-spiral pattern of Jordan's figures.

    P = zeros(nSteps+1, 2);
    P(1, :) = startPt;

    for t = 1:nSteps
        i_curr = mod(t-1, 4) + 1;
        i_next = mod(t,   4) + 1;

        C_curr = corners(i_curr, :);
        C_next = corners(i_next, :);

        delta     = P(t, :) - C_curr;        % current offset
        delta_rot = (R * delta')';           % rotate 90 deg CCW

        P(t+1, :) = C_next + gamma * delta_rot;
    end
end
