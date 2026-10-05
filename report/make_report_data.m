%% make_report_data.m -- regenerate every number and figure used in report.tex
%
% Writes:
%   report/figs/*.pdf        vector figures
%   report/data/*.tex        LaTeX fragments (matrices, tables) \input by report.tex
%
% Run from the report folder:  make_report_data

clear; clc; close all;

here  = fileparts(mfilename('fullpath'));
model = fullfile(here, '..');          % the model code is the repo root
figs  = fullfile(here, 'figs');
data  = fullfile(here, 'data');
if ~exist(figs,'dir'), mkdir(figs); end
if ~exist(data,'dir'), mkdir(data); end
addpath(model);

% Force light figures regardless of the MATLAB desktop theme. Without this a
% dark-theme session exports figures with a black background, which is unusable
% in a printed report.
set(groot, 'defaultFigureColor',   'w');
set(groot, 'defaultAxesColor',     'w');
set(groot, 'defaultAxesXColor',    'k');
set(groot, 'defaultAxesYColor',    'k');
set(groot, 'defaultAxesGridColor', [0.15 0.15 0.15]);
set(groot, 'defaultTextColor',     'k');
set(groot, 'defaultLegendColor',   'w');
set(groot, 'defaultLegendTextColor','k');
set(groot, 'defaultLegendEdgeColor','k');
cleanupTheme = onCleanup(@() reset(groot));

airframes = {'aerosonde','canard'};
R = struct();

%% ======================================================= the chain
for i = 1:numel(airframes)
    af = airframes{i};
    P  = lon_params(af);
    [xt, ut, info] = compute_lon_trim(P.Va0, 0, P);
    [A, B]  = compute_lon_ss(xt, ut, P);
    c       = compute_lon_tf(xt, ut, P);
    M       = lon_modes(A);
    AP      = lon_gains(c, P);

    R(i).af = af;  R(i).P = P;  R(i).xt = xt;  R(i).ut = ut;
    R(i).info = info;  R(i).A = A;  R(i).B = B;
    R(i).c = c;  R(i).M = M;  R(i).AP = AP;
end

%% ================================================ matrices as LaTeX
for i = 1:numel(R)
    writematrixtex(fullfile(data, sprintf('A_%s.tex', R(i).af)), R(i).A, '%.4f');
    writematrixtex(fullfile(data, sprintf('B_%s.tex', R(i).af)), R(i).B, '%.4f');
end

%% ============================================== trim + coef table
fid = fopen(fullfile(data,'trim_table.tex'),'w');
fprintf(fid, '\\begin{tabular}{lrr}\n\\toprule\n');
fprintf(fid, 'Quantity & Aerosonde & Canard \\\\\n\\midrule\n');
rows = {
 '$V_a^*$ [m/s]',          @(r) r.P.Va0,                  '%.2f'
 '$\alpha^*$ [deg]',       @(r) rad2deg(r.info.alpha),    '%.4f'
 '$\theta^*$ [deg]',       @(r) rad2deg(r.info.theta),    '%.4f'
 '$u^*$ [m/s]',            @(r) r.xt(1),                  '%.4f'
 '$w^*$ [m/s]',            @(r) r.xt(2),                  '%.4f'
 '$\delta_e^*$ [deg]',     @(r) rad2deg(r.ut(1)),         '%.4f'
 '$\delta_t^*$',           @(r) r.ut(2),                  '%.4f'
 '$\lVert f(x^*,u^*)\rVert$', @(r) r.info.residual,       '%.2e'
};
for k = 1:size(rows,1)
    fprintf(fid, '%s & $%s$ & $%s$ \\\\\n', rows{k,1}, ...
            sci2tex(num2str(rows{k,2}(R(1)), rows{k,3})), ...
            sci2tex(num2str(rows{k,2}(R(2)), rows{k,3})));
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);

fid = fopen(fullfile(data,'coef_table.tex'),'w');
fprintf(fid, '\\begin{tabular}{lrrr}\n\\toprule\n');
fprintf(fid, 'Coefficient & Aerosonde & Canard & Published ref. \\\\\n\\midrule\n');
ref = struct('a_theta1',5.294738359662443,'a_theta2',99.94742395724161, ...
             'a_theta3',-36.11239040790846,'a_V3',9.809999999999892);
cn  = {'a_theta1','a_theta2','a_theta3','a_V1','a_V2','a_V3'};
cl  = {'$a_{\theta_1}$','$a_{\theta_2}$','$a_{\theta_3}$', ...
       '$a_{V_1}$','$a_{V_2}$','$a_{V_3}$'};
for k = 1:numel(cn)
    if isfield(ref, cn{k})
        rs = sprintf('$%.6f$', ref.(cn{k}));
    else
        rs = '--';
    end
    fprintf(fid, '%s & $%.6f$ & $%.6f$ & %s \\\\\n', cl{k}, ...
            R(1).c.(cn{k}), R(2).c.(cn{k}), rs);
end
fprintf(fid, '\\midrule\n');
fprintf(fid, '$\\partial T/\\partial V_a$ & $%.4f$ & $%.4f$ & $-2.3522$ \\\\\n', ...
        R(1).c.dT_dVa, R(2).c.dT_dVa);
fprintf(fid, '$\\partial T/\\partial \\delta_t$ & $%.4f$ & $%.4f$ & $90.2794$ \\\\\n', ...
        R(1).c.dT_ddelta_t, R(2).c.dT_ddelta_t);
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);

%% ====================================================== mode table
fid = fopen(fullfile(data,'mode_table.tex'),'w');
fprintf(fid, '\\begin{tabular}{llrrr}\n\\toprule\n');
fprintf(fid, 'Airframe & Mode & $\\omega_n$ [rad/s] & $\\zeta$ & $T$ [s] \\\\\n\\midrule\n');
for i = 1:numel(R)
    nm = {'Short period (exact)','Short period (approx.)', ...
          'Phugoid (exact)','Phugoid (approx.)'};
    ss = {R(i).M.exact.short, R(i).M.approx.short, ...
          R(i).M.exact.phugoid, R(i).M.approx.phugoid};
    for k = 1:4
        if k == 1
            lbl = sprintf('\\texttt{%s}', R(i).af);
        else
            lbl = '';
        end
        fprintf(fid, '%s & %s & $%.4f$ & $%.4f$ & $%.3f$ \\\\\n', ...
                lbl, nm{k}, ss{k}.wn, ss{k}.zeta, ss{k}.period);
    end
    if i < numel(R), fprintf(fid, '\\midrule\n'); end
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);

%% ===================================================== gains table
fid = fopen(fullfile(data,'gain_table.tex'),'w');
fprintf(fid, '\\begin{tabular}{llrr}\n\\toprule\n');
fprintf(fid, 'Loop & Gain & Aerosonde & Canard \\\\\n\\midrule\n');
g = { 'Pitch (\S6.4.1)',    '$k_{p_\theta}$', 'kp_theta'
      '',                   '$k_{d_\theta}$', 'kd_theta'
      '',                   '$K_{\theta_{DC}}$','K_theta_DC'
      '',                   '$\omega_{n_\theta}$','wn_theta'
      'Altitude (\S6.4.2)', '$k_{p_h}$',      'kp_h'
      '',                   '$k_{i_h}$',      'ki_h'
      '$V_a$ on pitch (\S6.4.3)','$k_{p_{V_2}}$','kp_V2'
      '',                   '$k_{i_{V_2}}$',  'ki_V2'
      '$V_a$ on throttle (\S6.4.4)','$k_{p_V}$','kp_V'
      '',                   '$k_{i_V}$',      'ki_V' };
for k = 1:size(g,1)
    fprintf(fid, '%s & %s & $%.4f$ & $%.4f$ \\\\\n', g{k,1}, g{k,2}, ...
            R(1).AP.(g{k,3}), R(2).AP.(g{k,3}));
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);

%% ============================= fig 1: open-loop elevator doublet
P = R(1).P;  xt = R(1).xt;  ut = R(1).ut;
Ts = P.Ts;  N = round(30/Ts);  t = (0:N-1)'*Ts;
X = zeros(N,5);  U = zeros(N,2);
x = xt;
for k = 1:N
    de = ut(1);
    if t(k) >= 2 && t(k) < 2.5,       de = de + deg2rad(2);
    elseif t(k) >= 2.5 && t(k) < 3,   de = de - deg2rad(2);
    end
    uk = [de; ut(2)];
    X(k,:) = x';  U(k,:) = uk';
    x = rk4(x, uk, P, Ts);
end
dk = 1:5:N;   % decimate for plotting: vector export of 3000 pts/line is slow
f = figure('Position',[100 100 620 520]);
subplot(4,1,1); plot(t(dk), hypot(X(dk,1),X(dk,2)),'LineWidth',1.1); grid on;
ylabel('$V_a$ [m/s]','Interpreter','latex');
subplot(4,1,2); plot(t(dk), rad2deg(atan2(X(dk,2),X(dk,1))),'LineWidth',1.1); grid on;
ylabel('$\alpha$ [deg]','Interpreter','latex');
subplot(4,1,3); plot(t(dk), rad2deg(X(dk,4)),'LineWidth',1.1); grid on;
ylabel('$\theta$ [deg]','Interpreter','latex');
subplot(4,1,4); plot(t(dk), X(dk,5),'LineWidth',1.1); grid on;
ylabel('$h$ [m]','Interpreter','latex'); xlabel('$t$ [s]','Interpreter','latex');
exportgraphics(f, fullfile(figs,'openloop.pdf'), 'ContentType','vector', 'BackgroundColor','white');

%% ========================= fig 2: closed-loop altitude step
f = figure('Position',[100 100 620 560]);
for i = 1:numel(R)
    P = R(i).P;  xt = R(i).xt;  ut = R(i).ut;  AP = R(i).AP;
    Ts = P.Ts;  N = round(120/Ts);  t = (0:N-1)'*Ts;
    h_c = xt(5) + 50;  Va_c = P.Va0;
    X = zeros(N,5);  U = zeros(N,2);  TC = zeros(N,1);
    x = xt;  S = lon_autopilot_reset(ut(2));
    for k = 1:N
        [uk, S, dbg] = lon_autopilot([h_c; Va_c], x, S, AP);
        X(k,:) = x';  U(k,:) = uk';  TC(k) = dbg.theta_c;
        x = rk4(x, uk, P, Ts);
    end
    dk = 1:10:N;  % decimate for plotting
    subplot(4,2,i);   plot(t(dk),X(dk,5),'LineWidth',1.1); hold on;
    yline(h_c,'--'); grid on; title(sprintf('\\texttt{%s}',R(i).af),'Interpreter','latex');
    ylabel('$h$ [m]','Interpreter','latex');
    subplot(4,2,2+i); plot(t(dk),hypot(X(dk,1),X(dk,2)),'LineWidth',1.1); hold on;
    yline(Va_c,'--'); grid on; ylabel('$V_a$ [m/s]','Interpreter','latex');
    subplot(4,2,4+i); plot(t(dk),rad2deg(X(dk,4)),'LineWidth',1.1); hold on;
    plot(t(dk),rad2deg(TC(dk)),'--','LineWidth',0.8); grid on;
    ylabel('$\theta$ [deg]','Interpreter','latex');
    subplot(4,2,6+i); plot(t(dk),rad2deg(U(dk,1)),'LineWidth',1.1); hold on;
    plot(t(dk),U(dk,2)*10,'LineWidth',1.1); grid on;
    ylabel('cmd','Interpreter','latex'); xlabel('$t$ [s]','Interpreter','latex');
end
exportgraphics(f, fullfile(figs,'closedloop.pdf'), 'ContentType','vector', 'BackgroundColor','white');

%% ================================================ fig 3: pole map
% Two rows: the full plane, then a zoom on the phugoid. Without the zoom the
% phugoid pair sits within a pixel of the origin, because the short-period poles
% are roughly twenty times further out.
f = figure('Position',[100 100 660 480]);
cols = lines(2);
for i = 1:numel(R)
    lam = R(i).M.exact.all;
    la  = [R(i).M.approx.short.lambda; R(i).M.approx.phugoid.lambda];

    % --- full view
    subplot(2,2,i);
    plot(real(lam), imag(lam), 'o', 'MarkerSize',8, 'LineWidth',1.3, ...
         'Color', cols(1,:)); hold on;
    plot(real(la), imag(la), 'x', 'MarkerSize',10, 'LineWidth',1.3, ...
         'Color', cols(2,:));
    grid on; xline(0,'k-'); yline(0,'k-');
    xlabel('$\mathrm{Re}$','Interpreter','latex');
    ylabel('$\mathrm{Im}$','Interpreter','latex');
    title(sprintf('\\texttt{%s} --- full', R(i).af),'Interpreter','latex');
    if i == 1
        lg = legend({'$\mathrm{eig}(A_{\mathrm{lon}})$','\S5.6 approx.'}, ...
               'Interpreter','latex','Location','northwest','FontSize',7);
        set(lg, 'Color','w', 'TextColor','k', 'EdgeColor',[0.4 0.4 0.4]);
    end

    % --- phugoid zoom
    subplot(2,2,2+i);
    ph  = R(i).M.exact.phugoid.lambda;
    pha = R(i).M.approx.phugoid.lambda;
    plot(real(ph), imag(ph), 'o', 'MarkerSize',8, 'LineWidth',1.3, ...
         'Color', cols(1,:)); hold on;
    plot(real(pha), imag(pha), 'x', 'MarkerSize',10, 'LineWidth',1.3, ...
         'Color', cols(2,:));
    grid on; xline(0,'k-'); yline(0,'k-');
    span = max(abs([ph; pha]))*1.6;
    xlim([-span span]); ylim([-span span]);
    xlabel('$\mathrm{Re}$','Interpreter','latex');
    ylabel('$\mathrm{Im}$','Interpreter','latex');
    title(sprintf('\\texttt{%s} --- phugoid zoom', R(i).af),'Interpreter','latex');
end
exportgraphics(f, fullfile(figs,'poles.pdf'), 'ContentType','vector', ...
               'BackgroundColor','white');

%% ==================================================== console echo
fprintf('\nWrote figures to %s\n', figs);
fprintf('Wrote LaTeX data to %s\n', data);
for i = 1:numel(R)
    fprintf('  %-10s alpha*=%8.4f deg  de*=%8.4f deg  dt*=%.4f  |f|=%.2e\n', ...
            R(i).af, rad2deg(R(i).info.alpha), rad2deg(R(i).ut(1)), ...
            R(i).ut(2), R(i).info.residual);
end

%% ======================================================= helpers
function xn = rk4(x, u, P, dt)
k1 = f_lon(x,           u, P);
k2 = f_lon(x + dt/2*k1, u, P);
k3 = f_lon(x + dt/2*k2, u, P);
k4 = f_lon(x + dt*k3,   u, P);
xn = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

function s = sci2tex(s)
%SCI2TEX Turn MATLAB's 1.91e-15 into LaTeX 1.91\times10^{-15}.
% Needed because inside math mode a bare "e" renders as an italic variable.
s = regexprep(strtrim(s), '([0-9.]+)[eE]([+-]?)0*([0-9]+)', ...
              '$1\\times 10^{$2$3}');
end

function writematrixtex(fname, Mx, fmt)
%WRITEMATRIXTEX Write a matrix as a bare pmatrix body for \input.
fid = fopen(fname,'w');
for r = 1:size(Mx,1)
    for cc = 1:size(Mx,2)
        v = Mx(r,cc);
        if v == 0, s = '0'; else, s = sprintf(fmt, v); end
        if cc < size(Mx,2)
            fprintf(fid, '%s & ', s);
        else
            fprintf(fid, '%s', s);
        end
    end
    if r < size(Mx,1), fprintf(fid, ' \\\\\n'); else, fprintf(fid, '\n'); end
end
fclose(fid);
end
