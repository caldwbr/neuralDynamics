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

for Itest = 3:0.25:7
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

%% Quad window
fig = figure('Color','white', 'Position',[50 50 1100 900]);
plotOneRegime(p_SNIC, I_SNIC, subplot(2,2,1), 'SNIC');
plotOneRegime(p_SN,   I_SN,   subplot(2,2,2), 'SN');
plotOneRegime(p_subH, I_subH, subplot(2,2,3), 'Subcritical Hopf');
plotOneRegime(p_supH, I_supH, subplot(2,2,4), 'Supercritical Hopf');

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

function plotOneRegime(p, I, ax, regimeName)
    [dVdt, dndt, minf, ninf] = getModel(p, I);
    F = @(t,x) [dVdt(x(1),x(2)); dndt(x(1),x(2))];

    %% Plotting window (subHopf needs more n-headroom for its outer LC)
    Vmin = -90;  Vmax = 30;
    if strcmpi(regimeName,'Subcritical Hopf')
        nmin = -0.05;  nmax = 1.00;
    else
        nmin = -0.05;  nmax = 0.75;
    end

    %% Nullclines (USE minf for V-nullcline)
    Vfine = linspace(Vmin, Vmax, 1500);
    n_V = ( I - p.g_L*(Vfine - p.E_L) ...
             - p.g_Na*minf(Vfine).*(Vfine - p.E_Na) ) ...
          ./ ( p.g_K .* (Vfine - p.E_K) );          % V' = 0  =>  n = ...
    n_n = ninf(Vfine);                              % n' = 0
    mask = n_V > -0.1 & n_V < (nmax + 0.15);        % hide asymptotes

    %% Equilibria
    g_eq = @(V) I - p.g_L*(V - p.E_L) ...
                  - p.g_Na*minf(V).*(V - p.E_Na) ...
                  - p.g_K * ninf(V) .* (V - p.E_K);     % on n-nullcline
    diffN = n_V - n_n;
    idx   = find(diffN(1:end-1).*diffN(2:end) < 0);
    Veq = []; neq = [];
    for k = 1:numel(idx)
        try
            Vstar = fzero(g_eq, [Vfine(idx(k)) Vfine(idx(k)+1)]);
            Veq(end+1,1) = Vstar; neq(end+1,1) = ninf(Vstar);
        catch
        end
    end

    %% Classify equilibria via Jacobian
    h = 1e-4;
    stable_V=[]; stable_n=[]; unstable_V=[]; unstable_n=[]; saddle_V=[]; saddle_n=[];
    for k = 1:numel(Veq)
        V0 = Veq(k); n0 = neq(k);
        J = [(dVdt(V0+h,n0)-dVdt(V0-h,n0))/(2*h), (dVdt(V0,n0+h)-dVdt(V0,n0-h))/(2*h);
             (dndt(V0+h,n0)-dndt(V0-h,n0))/(2*h), (dndt(V0,n0+h)-dndt(V0,n0-h))/(2*h)];
        ev = eig(J); re = real(ev);
        if re(1)*re(2) < 0
            saddle_V(end+1,1)=V0; saddle_n(end+1,1)=n0;
        elseif all(re < 0)
            stable_V(end+1,1)=V0; stable_n(end+1,1)=n0;
        else
            unstable_V(end+1,1)=V0; unstable_n(end+1,1)=n0;
        end
    end

    axes(ax); hold on; box on;

    %% Streamlines (ode15s; NaN-break on box exit so plot doesn't connect across)
    opts = odeset('RelTol',1e-4,'AbsTol',1e-5);
    if p.tau < 0.2
        nSeedsV = 8; nSeedsN = 6; tmax = 15;
    else
        nSeedsV = 10; nSeedsN = 8; tmax = 25;
    end
    Vseeds = linspace(Vmin+10, Vmax-10, nSeedsV);
    Nseeds = linspace(nmin+0.05, nmax-0.10, nSeedsN);
    if ~isempty(saddle_V)
        th = linspace(0, 2*pi, 12);
        Vseeds = [Vseeds, saddle_V(1) + 6*cos(th)];
        Nseeds = [Nseeds, saddle_n(1) + 0.06*sin(th)];
    end
    for i = 1:length(Vseeds)
        for j = 1:length(Nseeds)
            if Vseeds(i)<Vmin||Vseeds(i)>Vmax||Nseeds(j)<nmin||Nseeds(j)>nmax, continue; end
            x0 = [Vseeds(i); Nseeds(j)];
            try
                [~,y] = ode15s(F, [0  tmax], x0, opts);
                oob = y(:,1)<Vmin | y(:,1)>Vmax | y(:,2)<nmin | y(:,2)>nmax;
                y(oob,:) = NaN;
                if sum(~oob) > 3
                    plot(ax,y(:,1),y(:,2),'Color',[0.55 0.55 0.6],'LineWidth',0.4,'HandleVisibility','off');
                end
            catch; end
            try
                [~,y] = ode15s(F, [0 -tmax], x0, opts);
                oob = y(:,1)<Vmin | y(:,1)>Vmax | y(:,2)<nmin | y(:,2)>nmax;
                y(oob,:) = NaN;
                if sum(~oob) > 3
                    plot(ax,y(:,1),y(:,2),'Color',[0.55 0.55 0.6],'LineWidth',0.4,'HandleVisibility','off');
                end
            catch; end
        end
    end

    %% Nullclines on top
    plot(ax, Vfine(mask), n_V(mask), 'b-', 'LineWidth',1.5, 'DisplayName','$\dot{V}=0$');
    plot(ax, Vfine, n_n,             'r-', 'LineWidth',1.5, 'DisplayName','$\dot{n}=0$');

    %% Limit cycle -- robust multi-seed (picks largest-amplitude closed orbit)
    opts_lc = odeset('RelTol',1e-5,'AbsTol',1e-5);
    seedList = [];
    if ~isempty(unstable_V)
        seedList = [seedList;  unstable_V(1)+1, unstable_n(1)+0.02];
    end
    if ~isempty(saddle_V)
        seedList = [seedList;  saddle_V(1)+2,   saddle_n(1)+0.05];
    end
    % Far-amplitude kicks -- catches outer LC of subHopf
    seedList = [seedList; ...
                -10, 0.55; ...
                  0, 0.40; ...
                -20, 0.65; ...
                -50, 0.05; ...
                 10, 0.85; ...      % outer subHopf LC
                -70, 0.10];         % left of focus, low n
    bestAmp = 0; X_best = [];
    for s = 1:size(seedList,1)
        x0s = seedList(s,:).';
        try
            [~,X1] = ode15s(F, [0 400], x0s, opts_lc);   % long settle
            [~,X2] = ode15s(F, [0  60], X1(end,:).', opts_lc);
            inBox = all(X2(:,1)>=Vmin-2 & X2(:,1)<=Vmax+2 & ...
                X2(:,2)>=nmin-0.05 & X2(:,2)<=nmax+0.05);
            if inBox && std(X2(:,1))>5
                amp = max(X2(:,1)) - min(X2(:,1));
                if amp > bestAmp
                    bestAmp = amp; X_best = X2;
                end
            end
        catch; end
    end
    if ~isempty(X_best)
        plot(ax,X_best(:,1),X_best(:,2),'Color',[0 0.6 0.2],'LineWidth',2, ...
             'DisplayName','limit cycle');
    end

    %% Separatrix (backward time, no events; clip after)
    if ~isempty(saddle_V)
        Vs = saddle_V(1); Ns = saddle_n(1);
        hJ = 1e-4;
        J = [(dVdt(Vs+hJ,Ns)-dVdt(Vs-hJ,Ns))/(2*hJ), (dVdt(Vs,Ns+hJ)-dVdt(Vs,Ns-hJ))/(2*hJ);
             (dndt(Vs+hJ,Ns)-dndt(Vs-hJ,Ns))/(2*hJ), (dndt(Vs,Ns+hJ)-dndt(Vs,Ns-hJ))/(2*hJ)];
        [Vc,Dc] = eig(J); ev = diag(Dc);
        [~,is] = min(real(ev));
        vs = Vc(:,is); vs = vs / norm(vs);
        eps0 = 5e-3;
        Fback = @(t,x) -F(t,x);
        opts_sep = odeset('RelTol',1e-5,'AbsTol',1e-5);
        for sgn = [+1 -1]
            try
                [~,S] = ode15s(Fback, [0 25], [Vs;Ns] + sgn*eps0*vs, opts_sep);
                oob = S(:,1)<Vmin|S(:,1)>Vmax|S(:,2)<nmin|S(:,2)>nmax;
                S(oob,:) = NaN;
                if sum(~oob) > 5
                    if sgn==+1
                        plot(ax,S(:,1),S(:,2),'Color',[0.55 0.1 0.55],'LineWidth',1.4,'LineStyle','--','DisplayName','separatrix');
                    else
                        plot(ax,S(:,1),S(:,2),'Color',[0.55 0.1 0.55],'LineWidth',1.4,'LineStyle','--','HandleVisibility','off');
                    end
                end
            catch; end
        end
    end

    %% Equilibria
    if ~isempty(stable_V)
        plot(ax,stable_V,stable_n,'ko','MarkerSize',10,'MarkerFaceColor',[0.23 0.55 0.94],'DisplayName','stable equilibrium');
    end
    if ~isempty(unstable_V)
        plot(ax,unstable_V,unstable_n,'ko','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor',[0.88 0.22 0.16],'LineWidth',1.5,'DisplayName','unstable equilibrium');
    end
    if ~isempty(saddle_V)
        plot(ax,saddle_V,saddle_n,'ks','MarkerSize',10,'MarkerFaceColor',[0.65 0.4 0.85],'DisplayName','saddle');
    end

    xlim(ax,[Vmin Vmax]); ylim(ax,[nmin nmax]);
    xlabel(ax,'$V$ (mV)','Interpreter','latex','FontSize',10);
    ylabel(ax,'$n$','Interpreter','latex','FontSize',10);
    title(ax,sprintf('%s, $I = %.1f$ pA', regimeName, I),'Interpreter','latex','FontSize',11);
    grid(ax,'on');
    legend(ax,'Location','best','Interpreter','latex','FontSize',7);
    hold(ax,'off');
end
