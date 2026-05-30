%% Classification of equilibria by trace tau and determinant Delta of L
%  Replicates the Izhikevich-style trace-determinant chart:
%    * Delta on horizontal axis, tau on vertical axis
%    * Parabola tau^2 - 4*Delta = 0
%    * Shaded region (Delta>0, tau<0) = asymptotic stability
%    * Six representative phase portraits with thick arrowheads
%    * Eigenvalue mini-diagrams beside each label
%
%  Phase portraits use a handful of full trajectories (4-8 each) with
%  filled-triangle arrowheads -- NOT vector-field quivers. Saddle shows
%  the stable + unstable manifolds explicitly along the eigenvectors,
%  plus a hyperbolic curve through each quadrant.

clear; clc; close all;

%% Plot ranges
Dmin = -2.2;  Dmax = 6.6;
Tmin = -5.0;  Tmax = 5.0;

%% Figure & main axes
fig = figure('Color','white','Position',[60 60 1400 900]);
axMainPos = [0.05 0.06 0.92 0.84];
ax = axes('Position', axMainPos); hold on;
xlim([Dmin Dmax]); ylim([Tmin Tmax]);

%% Shaded stable region: Delta > 0 AND tau < 0
patch([0 Dmax Dmax 0], [0 0 Tmin Tmin], [0.86 0.86 0.86], ...
      'EdgeColor','none','HandleVisibility','off');

%% Parabola Delta = tau^2 / 4  (i.e. tau^2 - 4*Delta = 0)
Tpar = linspace(Tmin, Tmax, 800);
Dpar = Tpar.^2 / 4;
plot(Dpar, Tpar, 'k-', 'LineWidth', 2.0);

%% Coordinate axes
plot([Dmin Dmax], [0 0], 'k-', 'LineWidth', 1.2);
plot([0 0], [Tmin Tmax], 'k-', 'LineWidth', 1.2);
% Arrow tips on axes
arrowTriTip(ax, [Dmax,        0], [1 0], 0.22, 0.10);
arrowTriTip(ax, [0,         Tmax], [0 1], 0.22, 0.10);

%% Axis-tip labels
text(Dmax-0.10, -0.35,    '\Delta', 'FontSize', 22, 'FontWeight','bold', 'FontAngle','italic', ...
     'HorizontalAlignment','right', 'VerticalAlignment','top');
%text(-0.15,     Tmax-0.15, '\tau',  'FontSize', 22, 'FontWeight','bold', 'FontAngle','italic', ...
%     'HorizontalAlignment','right', 'VerticalAlignment','top');
%text(-0.18, -0.32, '0', 'FontSize', 13);

%% Parabola equation labels (one on each branch)
idxTop = round(0.94*numel(Tpar));
idxBot = round(0.06*numel(Tpar));
text(Dpar(idxTop)+0.10, Tpar(idxTop)-0.25, '\tau^2 - 4\Delta = 0', 'FontSize', 26);
%text(Dpar(idxBot)+0.10, Tpar(idxBot)+0.25, '\tau^2 - 4\Delta = 0', 'FontSize', 26);

%% Bifurcation labels
text(0.15,  1.2, 'saddle-node bifurcation', 'FontSize', 22, ...
     'Rotation', 90, 'BackgroundColor','w');
text(0.15, -4.5, 'saddle-node bifurcation', 'FontSize', 22, ...
     'Rotation', 90, 'BackgroundColor','w');
text(2.8, 0.20, 'Andronov-Hopf bifurcation', 'FontSize', 24);

%% Hide MATLAB's default frame
set(ax, 'XColor','none','YColor','none','Color','none','Box','off');

%% Data -> normalized figure coordinate mapping
toFigX = @(D) axMainPos(1) + (D - Dmin)/(Dmax - Dmin) * axMainPos(3);
toFigY = @(T) axMainPos(2) + (T - Tmin)/(Tmax - Tmin) * axMainPos(4);
insetRect = @(D,T,w,h) [toFigX(D)-w/2, toFigY(T)-h/2, w, h];

W = 0.105; H = 0.140;   % phase-portrait inset size

%% --------- Regions: representative matrix + layout ---------
regions(1) = struct( ...                          % SADDLE
    'A',      [ 1.0, 0 ; 0, -1.0], ...
    'eigs',   [1.0, -1.0], ...
    'insetD', -1.05, 'insetT',  0.5, ...
    'labelD', -1.55, 'labelT', -0.7, ...
    'descD',  -1.55, 'descT',  -1.00, ...
    'eigD',   -1.15, 'eigT',   -1.80, ...
    'label',  'saddle', ...
    'desc',   '(real eigenvalues, different signs)', ...
    'type',   'saddle' );

regions(2) = struct( ...                          % UNSTABLE NODE
    'A',      [1.8, 0; 0, 1.0], ...
    'eigs',   [1.8, 1.0], ...
    'insetD', 0.9, 'insetT', 4.0, ...
    'labelD', 1.5, 'labelT', 3.65, ...
    'descD',  1.5, 'descT',  3.35, ...
    'eigD',   1.9, 'eigT',   4.40, ...
    'label',  'unstable node', ...
    'desc',   '(real positive eigenvalues)', ...
    'type',   'node' );

regions(3) = struct( ...                          % UNSTABLE FOCUS
    'A',      [0.4, -1.5; 1.5, 0.4], ...
    'eigs',   [0.4+1.5i, 0.4-1.5i], ...
    'insetD', 4.75, 'insetT',  1.85, ...
    'labelD', 5.40, 'labelT',  2.55, ...
    'descD',  5.40, 'descT',   2.15, ...
    'eigD',   5.80, 'eigT',    1.30, ...
    'label',  'unstable focus', ...
    'desc',   sprintf('(complex eigenvalues,\npositive real part)'), ...
    'type',   'focus' );

regions(4) = struct( ...                          % STABLE FOCUS
    'A',      [-0.4, -1.5; 1.5, -0.4], ...
    'eigs',   [-0.4+1.5i, -0.4-1.5i], ...
    'insetD', 4.75, 'insetT', -2.0, ...
    'labelD', 5.40, 'labelT', -1.25, ...
    'descD',  5.40, 'descT',  -1.65, ...
    'eigD',   5.80, 'eigT',   -2.50, ...
    'label',  'stable focus', ...
    'desc',   sprintf('(complex eigenvalues,\nnegative real part)'), ...
    'type',   'focus' );

regions(5) = struct( ...                          % STABLE NODE
    'A',      [-1.8, 0; 0, -1.0], ...
    'eigs',   [-1.8, -1.0], ...
    'insetD', 0.9, 'insetT', -4.0, ...
    'labelD', 1.5, 'labelT', -4.35, ...
    'descD',  1.5, 'descT',  -4.65, ...
    'eigD',   1.9, 'eigT',   -3.60, ...
    'label',  'stable node', ...
    'desc',   '(real negative eigenvalues)', ...
    'type',   'node' );

%% Labels + small eigenvalue diagrams (in main-axes coords)
for k = 1:numel(regions)
    r = regions(k);
    text(r.labelD, r.labelT, r.label, 'FontSize', 15, 'FontWeight','bold');
    if contains(r.desc, sprintf('\n'))
        text(r.descD, r.descT, strsplit(r.desc, sprintf('\n')), 'FontSize', 13);
    else
        text(r.descD, r.descT, r.desc, 'FontSize', 13);
    end
    drawEigDiagram(ax, r.eigD, r.eigT, 0.40, 0.55, r.eigs);
end

%% Phase-portrait insets
for k = 1:numel(regions)
    r = regions(k);
    pos = insetRect(r.insetD, r.insetT, W, H);
    axI = axes('Position', pos); hold(axI,'on'); box(axI,'on');
    drawPhasePortrait(axI, r.A, r.type, r.eigs);
    xlim(axI,[-2 2]); ylim(axI,[-2 2]);
    set(axI,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',1.0);
end

text(ax, -0.10, 4.5, char(964), 'FontName','Times New Roman', 'FontSize', 24, ...
     'FontWeight','bold', 'FontAngle','italic', 'Interpreter','none', ...
     'HorizontalAlignment','right', 'VerticalAlignment','middle', 'Clipping','off');

annotation(fig,'textbox',[0.05 0.91 0.92 0.06], ...
    'String','The \tau–\Delta atlas of planar equilibria', ...
    'FontSize',30,'FontWeight','bold','EdgeColor','none', ...
    'HorizontalAlignment','center');

%% ================== HELPERS ==================

function drawEigDiagram(ax, cD, cT, halfD, halfT, eigs)
    % Small Re/Im cross with filled dots at the eigenvalues.
    plot(ax, [cD-halfD, cD+halfD], [cT, cT], 'k-', 'LineWidth', 1.0);
    plot(ax, [cD, cD], [cT-halfT, cT+halfT], 'k-', 'LineWidth', 1.0);
    sxRe = 0.55 * halfD / 2.0;
    sxIm = 0.55 * halfT / 2.0;
    for j = 1:numel(eigs)
        edx = cD + real(eigs(j)) * sxRe;
        edy = cT + imag(eigs(j)) * sxIm;
        plot(ax, edx, edy, 'k.', 'MarkerSize', 18);
    end
end

function drawPhasePortrait(ax, A, type, eigs)
    % Faint coordinate cross
    plot(ax, [-2 2], [0 0], '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.5);
    plot(ax, [0 0], [-2 2], '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.5);
    switch type
        case 'saddle', drawSaddle(ax, A);
        case 'node',   drawNode(ax, A, all(real(eigs)<0));
        case 'focus',  drawFocus(ax, A, all(real(eigs)<0));
    end
    % Equilibrium marker (filled = stable, open = unstable)
    if all(real(eigs) < 0)
        plot(ax, 0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor','k');
    else
        plot(ax, 0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor','w', 'LineWidth', 1.0);
    end
end

function drawSaddle(ax, A)
    odefun = @(t,x) A*x;
    [V, D] = eig(A);
    evals = real(diag(D));
    [~, iU] = max(evals);   vU = V(:, iU); vU = vU/norm(vU);   % unstable
    [~, iS] = min(evals);   vS = V(:, iS); vS = vS/norm(vS);   % stable

    % Unstable manifold: 2 trajectories along +/- vU, flow OUT
    for sgn = [-1 1]
        x0 = sgn*vU*0.03;
        [~, Y] = ode45(odefun, [0 6], x0);
        Y = clipBox(Y, 1.85);
        if size(Y,1) > 2
            plot(ax, Y(:,1), Y(:,2), 'k-', 'LineWidth', 1.6);
            arrowOnPath(ax, Y, 0.80, 0.22);
        end
    end
    % Stable manifold: 2 trajectories along +/- vS, flow IN
    for sgn = [-1 1]
        x0 = sgn*vS*1.85;
        [~, Y] = ode45(odefun, [0 6], x0);
        Y = clipBox(Y, 1.85);
        if size(Y,1) > 2
            plot(ax, Y(:,1), Y(:,2), 'k-', 'LineWidth', 1.6);
            arrowOnPath(ax, Y, 0.25, 0.22);
        end
    end
    % Hyperbolic curves through each quadrant -- integrate both ways
    seeds = [1.2*vU + 1.2*vS, -1.2*vU + 1.2*vS, ...
            -1.2*vU - 1.2*vS,  1.2*vU - 1.2*vS]';
    for k = 1:size(seeds,1)
        x0 = seeds(k,:)';
        [~, Yf] = ode45(odefun, [0  5], x0);
        [~, Yb] = ode45(odefun, [0 -5], x0);
        Y = [flipud(Yb); Yf(2:end,:)];
        Y = clipBox(Y, 1.92);
        if size(Y,1) > 4
            plot(ax, Y(:,1), Y(:,2), 'k-', 'LineWidth', 1.3);
            arrowOnPath(ax, Y, 0.62, 0.20);
        end
    end
end

function drawNode(ax, A, isStable)
    odefun = @(t,x) A*x;
    angles = linspace(0, 2*pi, 9); angles(end) = [];   % 8 trajectories
    for ang = angles
        if isStable
            x0 = 1.85 * [cos(ang); sin(ang)];
            arrowFrac = 0.45;
        else
            x0 = 0.04 * [cos(ang); sin(ang)];
            arrowFrac = 0.75;
        end
        [~, Y] = ode45(odefun, [0 6], x0);
        Y = clipBox(Y, 1.9);
        if size(Y,1) > 2
            plot(ax, Y(:,1), Y(:,2), 'k-', 'LineWidth', 1.3);
            arrowOnPath(ax, Y, arrowFrac, 0.20);
        end
    end
end

function drawFocus(ax, A, isStable)
    odefun = @(t,x) A*x;
    if isStable
        starts = [1.75 0; 0 1.75; -1.75 0; 0 -1.75];
        tspan = [0 14];
        arrowFrac = 0.18;
    else
        starts = 0.05*[1 0; 0 1; -1 0; 0 -1];
        tspan = [0 9];
        arrowFrac = 0.92;
    end
    for k = 1:size(starts,1)
        x0 = starts(k,:)';
        [~, Y] = ode45(odefun, tspan, x0);
        Y = clipBox(Y, 1.9);
        if size(Y,1) > 5
            plot(ax, Y(:,1), Y(:,2), 'k-', 'LineWidth', 1.3);
            arrowOnPath(ax, Y, arrowFrac, 0.20);
        end
    end
end

function Y = clipBox(Y, lim)
    % Keep the first contiguous segment of Y that stays within [-lim, lim]^2.
    inBox = all(abs(Y) <= lim, 2);
    if isempty(inBox), return; end
    if ~inBox(1)
        i0 = find(inBox, 1, 'first');
        if isempty(i0), Y = []; return; end
        Y = Y(i0:end, :);
        inBox = all(abs(Y) <= lim, 2);
    end
    iEnd = find(~inBox, 1, 'first');
    if isempty(iEnd), iEnd = size(Y,1); else, iEnd = iEnd-1; end
    Y = Y(1:iEnd, :);
end

function arrowOnPath(ax, Y, frac, len)
    % Place a filled-triangle arrowhead at fractional position along Y.
    if size(Y,1) < 3, return; end
    mi = max(2, min(size(Y,1)-1, round(frac*size(Y,1))));
    d = Y(mi+1,:) - Y(mi-1,:);
    nrm = norm(d);
    if nrm < 1e-9, return; end
    d = d/nrm;
    tipPos = Y(mi,:) + d*len*0.35;     % put tip a bit forward of midpoint
    arrowTriTip(ax, tipPos, d, len, len*0.45);
end

function arrowTriTip(ax, tipPos, d, len, halfW)
    % Filled triangle with TIP at tipPos pointing in unit direction d.
    n = [-d(2), d(1)];
    base1 = tipPos - d*len + n*halfW;
    base2 = tipPos - d*len - n*halfW;
    patch(ax, [tipPos(1) base1(1) base2(1)], ...
              [tipPos(2) base1(2) base2(2)], 'k', ...
          'EdgeColor','k', 'LineWidth', 0.5);
end
