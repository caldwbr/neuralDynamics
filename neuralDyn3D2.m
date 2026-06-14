%% Four excitability regimes: SNIC, SN, subHopf, superHopf
%  Fixes for the subcritical-Hopf panel:
%    (1) Streamlines that leave the plotting box now break with NaN
%        instead of being deleted (no more horizontal connect-lines).
%    (2) Per-panel n-window: subHopf gets nmax=0.95 so the outer stable
%        LC fits inside the box (the LC reaches n ~ 0.85-0.9).
%    (3) Limit-cycle search gets a far-amplitude seed list and a longer
%        settling time so it actually lands on the outer LC.
%    (4) I_subH bumped slightly past the AH point for robustness.
clear; clc; close all;

%% Parameter sets for the four regimes
% SNIC (slow K+, saddle-node on invariant circle) -- Class 1 from Fig 4.1a
p_SNIC = struct('C',1, 'E_L',-80, 'E_Na',60, 'E_K',-90, ...
                'g_L',8, 'g_Na',20, 'g_K',10, ...
                'V_m_half',-20, 'k_m',15, ...
                'V_n_half',-25, 'k_n',5, 'tau',1.0);
I_SNIC = 2.0;

% SN (fast K+, saddle-node of equilibria; LC pre-exists)
p_SN   = struct('C',1, 'E_L',-80, 'E_Na',60, 'E_K',-90, ...
                'g_L',8, 'g_Na',20, 'g_K',10, ...
                'V_m_half',-20, 'k_m',15, ...
                'V_n_half',-25, 'k_n',5, 'tau',0.152);
I_SN   = 2.0;

for Itest = 2:0.25:7
    [dV,dn,minf,ninf] = getModel(p_SN, Itest);
    g_eq = @(V) Itest - p_SN.g_L*(V - p_SN.E_L) ...
              - p_SN.g_Na*minf(V).*(V - p_SN.E_Na) ...
              - p_SN.g_K*ninf(V).*(V - p_SN.E_K);
    Vs = linspace(-90,30,2000);
    crossings = sum(diff(sign(g_eq(Vs))) ~= 0);
    fprintf('I=%.2f: %d equilibria\n', Itest, crossings);
end

% Subcritical Hopf -- Fig 6.16 of Izhikevich (Fig 4.1b base + overrides)
p_subH = struct('C',1, 'E_L',-78, 'E_Na',60, 'E_K',-90, ...
                'g_L',1, 'g_Na',4, 'g_K',4, ...
                'V_m_half',-30, 'k_m',7, ...        % Fig 6.16 override
                'V_n_half',-45, 'k_n',5, 'tau',1.0);% Fig 4.1b values
I_subH = 40.0;    % nudged past the AH bifurcation for a robust outer LC

% Supercritical Hopf -- Class 2, Fig 4.1b base
p_supH = struct('C',1, 'E_L',-78, 'E_Na',60, 'E_K',-90, ...
                'g_L',8, 'g_Na',20, 'g_K',10, ...
                'V_m_half',-20, 'k_m',15, ...
                'V_n_half',-45, 'k_n',5, 'tau',1.0);
I_supH = 10.0;

probe(p_SNIC, 'SNIC',  2.0:0.25:7.0);
probe(p_SN,   'SN',    2.0:0.25:7.0);
probe(p_subH, 'subH', 40.0:1.00:55.0);
probe(p_supH, 'supH',  0.0:2.00:40.0);

%% =================== 3D BIFURCATION QUAD ===================
% Bifurcation locations from your diagnostic sweeps
Ibif_SNIC = 4.6;
Ibif_SN   = 4.6;
Ibif_subH = 48.9;
Ibif_supH = 14.6;

% Wider I-ranges spanning pre- and post-bifurcation
Ir_SNIC = [0.0   8.0];     % up to ~1.7x Ibif
Ir_SN   = [0.0   8.0];
Ir_subH = [35.0  60.0];
Ir_supH = [0.0   22.0];

% fig3D = figure('Color','white','Position',[60 60 1300 1000]);
% fprintf('\nRendering 3D panels (each takes ~15-30 s)...\n');
% plot3DRegime(p_SNIC, Ir_SNIC, Ibif_SNIC, subplot(2,2,1), 'SNIC');             fprintf('  SNIC done.\n');
% plot3DRegime(p_SN,   Ir_SN,   Ibif_SN,   subplot(2,2,2), 'SN');               fprintf('  SN done.\n');
% plot3DRegime(p_subH, Ir_subH, Ibif_subH, subplot(2,2,3), 'Subcritical Hopf'); fprintf('  subH done.\n');
% plot3DRegime(p_supH, Ir_supH, Ibif_supH, subplot(2,2,4), 'Supercritical Hopf'); fprintf('  supH done.\n');

% %% =================== LC-ONLY QUAD ===================
% figLC = figure('Color','white','Position',[80 80 1300 1000]);
% makeLCPanel(p_SNIC, [Ibif_SNIC,  Ir_SNIC(2)*1.2], Ibif_SNIC, subplot(2,2,1), 'SNIC LC manifold');
% makeLCPanel(p_SN,   [Ibif_SN*0.5, Ir_SN(2)*1.2],  Ibif_SN,   subplot(2,2,2), 'SN LC manifold');
% makeLCPanel(p_subH, [Ibif_subH*0.85, 65],         Ibif_subH, subplot(2,2,3), 'Subcritical Hopf LC manifold');
% makeLCPanel(p_supH, [Ibif_supH,  Ir_supH(2)*1.4], Ibif_supH, subplot(2,2,4), 'Supercritical Hopf LC manifold');
%% =================== LC + STRUCTURE QUAD ===================
figLC = figure('Color','white','Position',[80 80 1300 1000]);
makeLCPanel(p_SNIC, Ir_SNIC,   Ibif_SNIC, subplot(2,2,1), 'SNIC LC manifold');
makeLCPanel(p_SN,   Ir_SN,     Ibif_SN,   subplot(2,2,2), 'SN LC manifold');
makeLCPanel(p_subH, [35, 65],  Ibif_subH, subplot(2,2,3), 'Subcritical Hopf LC manifold');
makeLCPanel(p_supH, Ir_supH,   Ibif_supH, subplot(2,2,4), 'Supercritical Hopf LC manifold');

%% ============================ FUNCTIONS ============================
function [dVdt, dndt, minf, ninf] = getModel(p, I)
    minf = @(V) 1 ./ (1 + exp((p.V_m_half - V) / p.k_m));
    ninf = @(V) 1 ./ (1 + exp((p.V_n_half - V) / p.k_n));
    dVdt = @(V,n) (I - p.g_L*(V - p.E_L) ...
                   - p.g_Na*minf(V).*(V - p.E_Na) ...
                   - p.g_K*n.*(V - p.E_K)) / p.C;
    dndt = @(V,n) (ninf(V) - n) / p.tau;
end

function probe(p, name, Ilist)
    fprintf('\n%s sweep (look for neq jump or tr crossing zero):\n', name);
    fprintf('  %5s | neq | type    |   tr     det    | Re(lambda) +/- Im(lambda)\n','I');
    for Itest = Ilist
        [dV,dn,minf,ninf] = getModel(p, Itest);
        g_eq = @(V) Itest - p.g_L*(V - p.E_L) ...
                  - p.g_Na*minf(V).*(V - p.E_Na) ...
                  - p.g_K*ninf(V).*(V - p.E_K);
        Vs = linspace(-90,30,3000);
        idx = find(diff(sign(g_eq(Vs))) ~= 0);
        Veq = [];
        for k = 1:numel(idx)
            try
                Veq(end+1) = fzero(g_eq, [Vs(idx(k)) Vs(idx(k)+1)]); %#ok<AGROW>
            catch; end
        end
        nEq = numel(Veq);
        if nEq == 0
            fprintf('  %5.2f |  0  | -                                                \n', Itest);
            continue;
        end
        % Probe the highest-V equilibrium (the spiking-relevant one)
        V0 = max(Veq); n0 = ninf(V0);
        h = 1e-4;
        J = [(dV(V0+h,n0)-dV(V0-h,n0))/(2*h), (dV(V0,n0+h)-dV(V0,n0-h))/(2*h);
             (dn(V0+h,n0)-dn(V0-h,n0))/(2*h), (dn(V0,n0+h)-dn(V0,n0-h))/(2*h)];
        ev = eig(J); tr = trace(J); de = det(J);
        if imag(ev(1)) ~= 0
            kind = 'focus';
        elseif ev(1)*ev(2) < 0
            kind = 'saddle';
        else
            kind = 'node ';
        end
        fprintf('  %5.2f |  %d  | %s   | %+6.3f  %+6.3f | %+6.3f  +/- %5.3fi\n', ...
                Itest, nEq, kind, tr, de, real(ev(1)), abs(imag(ev(1))));
    end
end


function plot3DRegime(p, Irange, Ibif, ax, regimeName)
    axes(ax); hold on; grid on; box on; view(38, 22);
    Vmin = -90; Vmax = 30; nmin = -0.05;
    if strcmpi(regimeName,'Subcritical Hopf'), nmax = 1.00; else, nmax = 0.75; end

    %% Equilibrium branches
    Idense = linspace(Irange(1), Irange(2), 120);
    sB = []; uB = []; pB = [];
    for I = Idense
        [Ve, ne, knd] = findEqClass(p, I, Vmin, Vmax);
        for k = 1:numel(Ve)
            row = [Ve(k), ne(k), I];
            switch knd{k}
                case 'stable',   sB(end+1,:) = row;
                case 'unstable', uB(end+1,:) = row;
                case 'saddle',   pB(end+1,:) = row;
            end
        end
    end
    if ~isempty(sB), plot3(sB(:,1),sB(:,2),sB(:,3),'-', 'Color',[0.05 0.20 0.55],'LineWidth',2.8,'DisplayName','stable equilibrium branch'); end
    if ~isempty(uB), plot3(uB(:,1),uB(:,2),uB(:,3),':', 'Color',[0.70 0.05 0.10],'LineWidth',2.8,'DisplayName','unstable equilibrium branch'); end
    if ~isempty(pB), plot3(pB(:,1),pB(:,2),pB(:,3),'--','Color',[0.45 0.15 0.60],'LineWidth',2.4,'DisplayName','saddle branch'); end

    %% 30 slices with nullclines + basin dots (only first slice contributes to legend)
    Islices = linspace(Irange(1), Irange(2), 30);
    for k = 1:numel(Islices)
        addShadedSlice3D(p, Islices(k), Vmin, Vmax, nmin, nmax, k==1);
    end

    %% LC manifold (stacked rings) — uses near-focus seeds for small-amp Hopf cycles
    plotLCManifold(p, linspace(Irange(1), Irange(2), 50), Vmin, Vmax, nmin, nmax);

    %% Bifurcation marker
    Iprobe = max(Ibif - 0.05, Irange(1) + 0.01);
    [Vep, nep, knp] = findEqClass(p, Iprobe, Vmin, Vmax);
    stI = find(strcmp(knp,'stable')); sdI = find(strcmp(knp,'saddle'));
    if ~isempty(stI) && ~isempty(sdI)
        Vmark = (Vep(stI(1)) + Vep(sdI(1)))/2;
        nmark = (nep(stI(1)) + nep(sdI(1)))/2;
    else
        [Vex, nex] = findEqClass(p, Ibif, Vmin, Vmax);
        Vmark = mean(Vex); nmark = mean(nex);
    end
    plot3(Vmark, nmark, Ibif, 'p', 'MarkerSize', 28, ...
          'MarkerFaceColor',[1 0.85 0.1],'MarkerEdgeColor','k','LineWidth',1.5, ...
          'DisplayName', sprintf('bifurcation, I=%.2f', Ibif));

    xlim([Vmin Vmax]); ylim([nmin nmax]); zlim(Irange);
    xlabel('V (mV)','FontSize',10); ylabel('n','FontSize',10); zlabel('I (pA)','FontSize',10);
    title(sprintf('%s (I_{bif} \\approx %.2f)', regimeName, Ibif),'FontSize',11);
    legend('Location','northeastoutside','FontSize',8);
    set(ax,'Projection','perspective');
end

function addShadedSlice3D(p, I, Vmin, Vmax, nmin, nmax, firstSlice)
    [dV,dn,minf,ninf] = getModel(p, I);
    F = @(t,x) [dV(x(1),x(2)); dn(x(1),x(2))];

    %% Nullclines at z = I
    Vfine = linspace(Vmin, Vmax, 600);
    n_V = (I - p.g_L*(Vfine - p.E_L) ...
             - p.g_Na*minf(Vfine).*(Vfine - p.E_Na)) ./ (p.g_K.*(Vfine - p.E_K));
    mask = n_V > nmin-0.1 & n_V < nmax+0.1;
    vCol = [0.15 0.55 0.65];  % teal
    nCol = [0.95 0.55 0.15];  % orange
    if firstSlice
        plot3(Vfine(mask), n_V(mask),  I*ones(sum(mask),1), '-', 'Color', vCol, 'LineWidth', 0.7, 'DisplayName','V-nullcline ($\dot V=0$)');
        plot3(Vfine,       ninf(Vfine), I*ones(numel(Vfine),1), '-', 'Color', nCol, 'LineWidth', 0.7, 'DisplayName','n-nullcline ($\dot n=0$)');
    else
        plot3(Vfine(mask), n_V(mask),  I*ones(sum(mask),1), '-', 'Color', vCol, 'LineWidth', 0.7, 'HandleVisibility','off');
        plot3(Vfine,       ninf(Vfine), I*ones(numel(Vfine),1), '-', 'Color', nCol, 'LineWidth', 0.7, 'HandleVisibility','off');
    end

    %% Basin sampling (small dots)
    [Ve, ne, knd] = findEqClass(p, I, Vmin, Vmax);
    stIdx = find(strcmp(knd,'stable'));
    nGV = 36; nGN = 22;
    Vg = linspace(Vmin+1, Vmax-1, nGV);
    Ng = linspace(nmin+0.01, nmax-0.01, nGN);
    [VV, NN] = meshgrid(Vg, Ng);
    C = ones(size(VV));
    if ~isempty(stIdx)
        opts = odeset('RelTol',1e-2,'AbsTol',1e-3);
        for ii = 1:numel(VV)
            try
                [~,y] = ode23(F, [0 30], [VV(ii); NN(ii)], opts);
                ep = y(end,:);
                for s = stIdx
                    if abs(ep(1)-Ve(s)) < 5 && abs(ep(2)-ne(s)) < 0.06
                        C(ii) = 0; break;
                    end
                end
            catch
            end
        end
    end
    cIn  = [0.55 0.75 0.60];   % sage
    cOut = [0.85 0.65 0.60];   % dusty rose
    inMask  = C == 0;
    outMask = C == 1;
    if any(inMask(:))
        if firstSlice
            plot3(VV(inMask), NN(inMask), I*ones(sum(inMask(:)),1), '.', 'Color', cIn,  'MarkerSize', 3.5, 'DisplayName','rest-basin sample');
        else
            plot3(VV(inMask), NN(inMask), I*ones(sum(inMask(:)),1), '.', 'Color', cIn,  'MarkerSize', 3.5, 'HandleVisibility','off');
        end
    end
    if any(outMask(:))
        if firstSlice
            plot3(VV(outMask), NN(outMask), I*ones(sum(outMask(:)),1), '.', 'Color', cOut, 'MarkerSize', 3.5, 'DisplayName','escape/LC-basin sample');
        else
            plot3(VV(outMask), NN(outMask), I*ones(sum(outMask(:)),1), '.', 'Color', cOut, 'MarkerSize', 3.5, 'HandleVisibility','off');
        end
    end
end

function [Veq, neq, kinds] = findEqClass(p, I, Vmin, Vmax)
    [dV,dn,minf,ninf] = getModel(p, I);
    g = @(V) I - p.g_L*(V - p.E_L) ...
              - p.g_Na*minf(V).*(V - p.E_Na) ...
              - p.g_K*ninf(V).*(V - p.E_K);
    Vf = linspace(Vmin, Vmax, 2000);
    idx = find(diff(sign(g(Vf))) ~= 0);
    Veq = []; neq = []; kinds = {}; h = 1e-4;
    for k = 1:numel(idx)
        try
            Vs = fzero(g, [Vf(idx(k)) Vf(idx(k)+1)]);
            ns = ninf(Vs);
            J = [(dV(Vs+h,ns)-dV(Vs-h,ns))/(2*h), (dV(Vs,ns+h)-dV(Vs,ns-h))/(2*h);
                 (dn(Vs+h,ns)-dn(Vs-h,ns))/(2*h), (dn(Vs,ns+h)-dn(Vs,ns-h))/(2*h)];
            ev = eig(J); re = real(ev);
            if re(1)*re(2) < 0,    knd = 'saddle';
            elseif all(re < 0),    knd = 'stable';
            else,                  knd = 'unstable';
            end
            Veq(end+1) = Vs; neq(end+1) = ns; kinds{end+1} = knd;
        catch; end
    end
end

function plotLCManifold(p, Ilist, Vmin, Vmax, nmin, nmax)
    baseSeeds = [10 0.5; -10 0.6; 0 0.4; -20 0.7; 10 0.85; -50 0.05; -70 0.10];
    opts = odeset('RelTol',1e-6,'AbsTol',1e-7);
    lcCol = [0.90 0.15 0.60];
    legShown = false;
    for I = Ilist
        [dV,dn,~,~] = getModel(p, I);
        F = @(t,x) [dV(x(1),x(2)); dn(x(1),x(2))];

        % Augment seeds with near-unstable-focus kicks (catches small-amp Hopf LCs)
        [Ve, ne, knd] = findEqClass(p, I, Vmin, Vmax);
        unstIdx = find(strcmp(knd,'unstable'));
        seeds = baseSeeds;
        for u = unstIdx
            seeds = [seeds; Ve(u)+0.2, ne(u)+0.003;
                            Ve(u)+0.8, ne(u)+0.015;
                            Ve(u)+2.0, ne(u)+0.04];
        end

        bestAmp = 0; LC = [];
        for s = 1:size(seeds,1)
            try
                [~,X1] = ode15s(F, [0 600], seeds(s,:).', opts);
                [~,X2] = ode15s(F, [0 120], X1(end,:).', opts);
                inBox = all(X2(:,1)>=Vmin-2 & X2(:,1)<=Vmax+2 & ...
                            X2(:,2)>=nmin-0.05 & X2(:,2)<=nmax+0.05);
                if inBox && std(X2(:,1)) > 1.5   % loose enough to catch small Hopf cycles
                    amp = max(X2(:,1)) - min(X2(:,1));
                    if amp > bestAmp, bestAmp = amp; LC = X2; end
                end
            catch; end
        end
        if ~isempty(LC)
            if ~legShown
                plot3(LC(:,1), LC(:,2), I*ones(size(LC,1),1), '-', ...
                      'Color',lcCol,'LineWidth',1.8,'DisplayName','limit cycle');
                legShown = true;
            else
                plot3(LC(:,1), LC(:,2), I*ones(size(LC,1),1), '-', ...
                      'Color',lcCol,'LineWidth',1.8,'HandleVisibility','off');
            end
        end
    end
end

function makeLCPanel(p, Irange, Ibif, ax, name)
    axes(ax); hold on; grid on; box on; view(40, 22);
    Vmin = -90; Vmax = 30; nmin = -0.05;
    if contains(lower(name),'subcritical'), nmax = 1.00; else, nmax = 0.75; end

    %% Equilibrium branches (same colors as bifurcation quad)
    Idense = linspace(Irange(1), Irange(2), 100);
    sB = []; uB = []; pB = [];
    for I = Idense
        [Ve, ne, knd] = findEqClass(p, I, Vmin, Vmax);
        for k = 1:numel(Ve)
            row = [Ve(k), ne(k), I];
            switch knd{k}
                case 'stable',   sB(end+1,:) = row; %#ok<*AGROW>
                case 'unstable', uB(end+1,:) = row;
                case 'saddle',   pB(end+1,:) = row;
            end
        end
    end
    if ~isempty(sB), plot3(sB(:,1),sB(:,2),sB(:,3),'-', 'Color',[0.05 0.20 0.55],'LineWidth',2.6,'DisplayName','stable equilibrium branch'); end
    if ~isempty(uB), plot3(uB(:,1),uB(:,2),uB(:,3),':', 'Color',[0.70 0.05 0.10],'LineWidth',2.6,'DisplayName','unstable equilibrium branch'); end
    if ~isempty(pB), plot3(pB(:,1),pB(:,2),pB(:,3),'--','Color',[0.45 0.15 0.60],'LineWidth',2.2,'DisplayName','saddle branch'); end

    %% Separatrix manifold (deep purple) — wherever a saddle exists
    Ilist_sep = linspace(Irange(1), Irange(2), 35);
    [XL,YL,ZL, XR,YR,ZR] = buildSeparatrixSurface(p, Ilist_sep, Vmin, Vmax, nmin, nmax);
    sepCol = [0.30 0.10 0.45];
    legSepShown = false;
    goodL = ~all(isnan(XL),2);
    if sum(goodL) >= 2
        surf(XL(goodL,:), YL(goodL,:), ZL(goodL,:), ...
             'FaceColor', sepCol, 'EdgeColor','none', 'FaceAlpha', 0.40, ...
             'FaceLighting','gouraud', 'DisplayName','separatrix manifold');
        legSepShown = true;
    end
    goodR = ~all(isnan(XR),2);
    if sum(goodR) >= 2
        nameArg = {'HandleVisibility','off'};
        if ~legSepShown, nameArg = {'DisplayName','separatrix manifold'}; end
        surf(XR(goodR,:), YR(goodR,:), ZR(goodR,:), ...
             'FaceColor', sepCol, 'EdgeColor','none', 'FaceAlpha', 0.40, ...
             'FaceLighting','gouraud', nameArg{:});
    end

    %% Unstable (saddle) limit-cycle manifold — present in subcritical-Hopf bistable region
    Ilist_ulc = linspace(Irange(1), Irange(2), 35);
    [Xu, Yu, Zu] = buildUnstableLCSurface(p, Ilist_ulc, Vmin, Vmax, nmin, nmax);
    goodu = ~all(isnan(Xu),2);
    if sum(goodu) >= 2
        surf(Xu(goodu,:), Yu(goodu,:), Zu(goodu,:), ...
             'FaceColor',[0.70 0.05 0.10], 'EdgeColor','none', ...
             'FaceAlpha', 0.45, 'FaceLighting','gouraud', ...
             'DisplayName','unstable (saddle) limit cycle');
    end

    %% LC manifold (magenta)
    Ilist_lc = linspace(Irange(1), Irange(2), 45);
    [X, Y, Z] = buildLCSurface(p, Ilist_lc, Vmin, Vmax, nmin, nmax);
    good = ~all(isnan(X),2);
    if sum(good) >= 2
        surf(X(good,:), Y(good,:), Z(good,:), 'FaceColor',[0.90 0.15 0.60], ...
             'EdgeColor','none','FaceAlpha',0.55,'FaceLighting','gouraud', ...
             'DisplayName','limit-cycle manifold');
    end

    if (sum(goodL)>=2 || sum(goodR)>=2 || sum(good)>=2)
        camlight headlight; lighting gouraud;
    end

    % Bifurcation level plane
    patch('XData',[Vmin Vmax Vmax Vmin],'YData',[nmin nmin nmax nmax], ...
          'ZData',Ibif*[1 1 1 1],'FaceColor',[1 0.85 0.1], ...
          'FaceAlpha',0.10,'EdgeColor',[0.6 0.5 0],'LineStyle',':', ...
          'HandleVisibility','off');

    xlim([Vmin Vmax]); ylim([nmin nmax]); zlim(Irange);
    xlabel('V (mV)'); ylabel('n'); zlabel('I (pA)');
    title(name,'FontSize',11);
    legend('Location','northeast','FontSize',8);
    set(ax,'Projection','perspective');
end

function [XL, YL, ZL, XR, YR, ZR] = buildSeparatrixSurface(p, Ilist, Vmin, Vmax, nmin, nmax)
    nArc = 60;
    XL = nan(numel(Ilist), nArc); YL = nan(numel(Ilist), nArc); ZL = nan(numel(Ilist), nArc);
    XR = nan(numel(Ilist), nArc); YR = nan(numel(Ilist), nArc); ZR = nan(numel(Ilist), nArc);
    opts = odeset('RelTol',1e-6,'AbsTol',1e-7);

    for kk = 1:numel(Ilist)
        I = Ilist(kk);
        [Ve, ne, knd] = findEqClass(p, I, Vmin, Vmax);
        sIdx = find(strcmp(knd,'saddle'));
        if isempty(sIdx), continue; end

        [dV,dn,~,~] = getModel(p, I);
        F = @(t,x) [dV(x(1),x(2)); dn(x(1),x(2))];

        Vs = Ve(sIdx(1)); Ns = ne(sIdx(1));
        h = 1e-4;
        J = [(dV(Vs+h,Ns)-dV(Vs-h,Ns))/(2*h), (dV(Vs,Ns+h)-dV(Vs,Ns-h))/(2*h);
             (dn(Vs+h,Ns)-dn(Vs-h,Ns))/(2*h), (dn(Vs,Ns+h)-dn(Vs,Ns-h))/(2*h)];
        [Vc, Dc] = eig(J); ev = diag(Dc);
        [~, is] = min(real(ev));
        vS = Vc(:,is)/norm(Vc(:,is));
        Fb = @(t,x) -F(t,x);
        eps0 = 5e-3;

        for sgn = [+1 -1]
            try
                [~, S] = ode15s(Fb, [0 25], [Vs;Ns] + sgn*eps0*vS, opts);
                oob = S(:,1)<Vmin | S(:,1)>Vmax | S(:,2)<nmin | S(:,2)>nmax;
                firstOob = find(oob, 1);
                if ~isempty(firstOob), S = S(1:firstOob-1,:); end
                if size(S,1) > 5
                    arcLen = [0; cumsum(sqrt(diff(S(:,1)).^2 + diff(S(:,2)).^2))];
                    arcQ = linspace(0, arcLen(end), nArc);
                    Xq = interp1(arcLen, S(:,1), arcQ);
                    Yq = interp1(arcLen, S(:,2), arcQ);
                    if sgn == +1
                        XL(kk,:) = Xq; YL(kk,:) = Yq; ZL(kk,:) = I;
                    else
                        XR(kk,:) = Xq; YR(kk,:) = Yq; ZR(kk,:) = I;
                    end
                end
            catch; end
        end
    end
end

function [X, Y, Z] = buildLCSurface(p, Ilist, Vmin, Vmax, nmin, nmax)
    nPhase = 80;
    X = nan(numel(Ilist), nPhase);
    Y = nan(numel(Ilist), nPhase);
    Z = nan(numel(Ilist), nPhase);

    baseSeeds = [10 0.5; -10 0.6; 0 0.4; -20 0.7; 10 0.85; -50 0.05; -70 0.10];
    opts = odeset('RelTol',1e-7,'AbsTol',1e-8);

    for kk = 1:numel(Ilist)
        I = Ilist(kk);
        [dV,dn,~,~] = getModel(p, I);
        F = @(t,x) [dV(x(1),x(2)); dn(x(1),x(2))];

        [Ve, ne, knd] = findEqClass(p, I, Vmin, Vmax);
        unstIdx = find(strcmp(knd,'unstable'));
        seeds = baseSeeds;
        for u = unstIdx
            seeds = [seeds; Ve(u)+0.2, ne(u)+0.003;
                            Ve(u)+0.8, ne(u)+0.015;
                            Ve(u)+2.0, ne(u)+0.04]; %#ok<*AGROW>
        end

        bestAmp = 0; LC = [];
        for s = 1:size(seeds,1)
            try
                [~,X1] = ode15s(F, [0 600], seeds(s,:).', opts);
                [~,X2] = ode15s(F, [0 120], X1(end,:).', opts);
                inBox = all(X2(:,1)>=Vmin-2 & X2(:,1)<=Vmax+2 & ...
                            X2(:,2)>=nmin-0.05 & X2(:,2)<=nmax+0.05);
                if inBox && std(X2(:,1)) > 0.8
                    amp = max(X2(:,1)) - min(X2(:,1));
                    if amp > bestAmp, bestAmp = amp; LC = X2; end
                end
            catch; end
        end

        if ~isempty(LC) && size(LC,1) > 20
            cx = mean(LC(:,1)); cy = mean(LC(:,2));
            sx = max(LC(:,1)) - min(LC(:,1)) + eps;
            sy = max(LC(:,2)) - min(LC(:,2)) + eps;
            ang = atan2((LC(:,2)-cy)/sy, (LC(:,1)-cx)/sx);
            [angSort, ix] = sort(ang);
            LCs = LC(ix,:);
            [angU, uix] = unique(angSort);
            if numel(angU) > 20
                LCu = LCs(uix,:);
                phaseQ = linspace(-pi, pi, nPhase);
                X(kk,:) = interp1(angU, LCu(:,1), phaseQ, 'linear', 'extrap');
                Y(kk,:) = interp1(angU, LCu(:,2), phaseQ, 'linear', 'extrap');
                Z(kk,:) = I;
            end
        end
    end
end

function [X, Y, Z] = buildUnstableLCSurface(p, Ilist, Vmin, Vmax, nmin, nmax)
    nPhase = 80;
    X = nan(numel(Ilist), nPhase);
    Y = nan(numel(Ilist), nPhase);
    Z = nan(numel(Ilist), nPhase);
    opts = odeset('RelTol',1e-7,'AbsTol',1e-8);

    for kk = 1:numel(Ilist)
        I = Ilist(kk);
        [Ve, ne, knd] = findEqClass(p, I, Vmin, Vmax);
        stIdx = find(strcmp(knd,'stable'));
        if isempty(stIdx), continue; end

        [dV,dn,~,~] = getModel(p, I);
        F  = @(t,x) [dV(x(1),x(2)); dn(x(1),x(2))];
        Fb = @(t,x) -F(t,x);

        Vstar = Ve(stIdx(1)); nstar = ne(stIdx(1));
        % Seeds slightly off the stable focus — backward-time spirals outward to the saddle cycle
        deltas = [0.3 0.003; 0.8 0.008; 1.5 0.020; 3.0 0.045; 5.0 0.080];

        bestAmp = 0; ULC = [];
        for d = 1:size(deltas,1)
            try
                x0 = [Vstar + deltas(d,1); nstar + deltas(d,2)];
                [~,X1] = ode15s(Fb, [0 400], x0, opts);
                [~,X2] = ode15s(Fb, [0 80], X1(end,:).', opts);
                inBox = all(X2(:,1) >= Vmin-2 & X2(:,1) <= Vmax+2 & ...
                            X2(:,2) >= nmin-0.05 & X2(:,2) <= nmax+0.05);
                if inBox && std(X2(:,1)) > 1.0
                    amp = max(X2(:,1)) - min(X2(:,1));
                    if amp > bestAmp, bestAmp = amp; ULC = X2; end
                end
            catch; end
        end

        if ~isempty(ULC) && size(ULC,1) > 20
            cx = mean(ULC(:,1)); cy = mean(ULC(:,2));
            sx = max(ULC(:,1)) - min(ULC(:,1)) + eps;
            sy = max(ULC(:,2)) - min(ULC(:,2)) + eps;
            ang = atan2((ULC(:,2)-cy)/sy, (ULC(:,1)-cx)/sx);
            [angSort, ix] = sort(ang);
            ULCs = ULC(ix,:);
            [angU, uix] = unique(angSort);
            if numel(angU) > 20
                ULCu = ULCs(uix,:);
                phaseQ = linspace(-pi, pi, nPhase);
                X(kk,:) = interp1(angU, ULCu(:,1), phaseQ, 'linear', 'extrap');
                Y(kk,:) = interp1(angU, ULCu(:,2), phaseQ, 'linear', 'extrap');
                Z(kk,:) = I;
            end
        end
    end
end
