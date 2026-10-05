%% run_lon.m -- longitudinal design chain, Beard & McLain chapters 3 to 6
%
% Edit `airframe` below to switch aircraft. Everything downstream adapts.
%
%   ch 3-4  nonlinear 5-state longitudinal plant      f_lon
%   ch 5    trim, linearisation, transfer functions, modes
%   ch 6    autopilot gains and a closed-loop altitude step
%
% Run:  run_lon

clear; clc; close all;

airframe = 'aerosonde';          % 'aerosonde' or 'canard'
Va_star  = 25;                   % m/s
gamma_star = 0;                  % rad

here = fileparts(mfilename('fullpath'));
addpath(here);

%% ===================================================== ch 3-4: the plant
P = lon_params(airframe);
fprintf('\n===== %s, Va* = %g m/s, gamma* = %g deg =====\n', ...
        upper(airframe), Va_star, rad2deg(gamma_star));

%% ========================================================= ch 5: trim
[x_trim, u_trim, info] = compute_lon_trim(Va_star, gamma_star, P);

fprintf('\n--- Trim (book sec 5.3 / appendix F.2) ---\n');
fprintf('  alpha*   = %10.6f rad  (%7.3f deg)\n', info.alpha, rad2deg(info.alpha));
fprintf('  theta*   = %10.6f rad  (%7.3f deg)\n', info.theta, rad2deg(info.theta));
fprintf('  u*, w*   = %10.6f, %9.6f m/s\n', x_trim(1), x_trim(2));
fprintf('  delta_e* = %10.6f rad  (%7.3f deg)\n', u_trim(1), rad2deg(u_trim(1)));
fprintf('  delta_t* = %10.6f\n', u_trim(2));
fprintf('  residual |f(x*,u*)| = %.3e   climb rate = %.4f m/s\n', ...
        info.residual, info.hdot);

%% ============================================ ch 5: linear state space
[A_lon, B_lon] = compute_lon_ss(x_trim, u_trim, P);

fprintf('\n--- A_lon  [u w q theta h]  (book eq 5.50) ---\n');
disp(A_lon);
fprintf('--- B_lon  [delta_e delta_t] ---\n');
disp(B_lon);

%% ========================================== ch 5: transfer functions
coef = compute_lon_tf(x_trim, u_trim, P);

fprintf('--- Transfer function coefficients (book sec 5.4) ---\n');
fprintf('  a_theta1 = %12.6f   a_V1 = %12.6f\n', coef.a_theta1, coef.a_V1);
fprintf('  a_theta2 = %12.6f   a_V2 = %12.6f\n', coef.a_theta2, coef.a_V2);
fprintf('  a_theta3 = %12.6f   a_V3 = %12.6f\n', coef.a_theta3, coef.a_V3);
fprintf('  dT/dVa   = %12.6f   dT/ddelta_t = %9.4f\n', ...
        coef.dT_dVa, coef.dT_ddelta_t);

fprintf('\n');
show_tf('theta / delta_e', coef.T_theta_delta_e);   % eq 5.29
show_tf('h / theta      ', coef.T_h_theta);         % eq 5.31
show_tf('h / Va         ', coef.T_h_Va);            % eq 5.32
show_tf('Va / delta_t   ', coef.T_Va_delta_t);      % eq 5.36
show_tf('Va / theta     ', coef.T_Va_theta);        % eq 5.36

%% ==================================================== ch 5: modes
M = lon_modes(A_lon);

fprintf('\n--- Reduced-order modes (book sec 5.6) ---\n');
fprintf('  %-14s %10s %10s %10s\n', '', 'wn rad/s', 'zeta', 'period s');
fprintf('  %-14s %10.4f %10.4f %10.3f\n', 'short exact', ...
        M.exact.short.wn,   M.exact.short.zeta,   M.exact.short.period);
fprintf('  %-14s %10.4f %10.4f %10.3f\n', 'short approx', ...
        M.approx.short.wn,  M.approx.short.zeta,  M.approx.short.period);
fprintf('  %-14s %10.4f %10.4f %10.3f\n', 'phugoid exact', ...
        M.exact.phugoid.wn, M.exact.phugoid.zeta, M.exact.phugoid.period);
fprintf('  %-14s %10.4f %10.4f %10.3f\n', 'phugoid approx', ...
        M.approx.phugoid.wn, M.approx.phugoid.zeta, M.approx.phugoid.period);

%% ============================================ persist for chapter 6
save(fullfile(here,'lon_coef.mat'), ...
     'coef','A_lon','B_lon','x_trim','u_trim','airframe');

%% ==================================================== ch 6: gains
AP = lon_gains(coef, P);

fprintf('\n--- Autopilot gains (book sec 6.4) ---\n');
fprintf('  pitch     kp = %9.4f  kd = %9.4f  K_DC = %7.4f  wn = %6.3f\n', ...
        AP.kp_theta, AP.kd_theta, AP.K_theta_DC, AP.wn_theta);
fprintf('  altitude  kp = %9.4f  ki = %9.4f               wn = %6.3f\n', ...
        AP.kp_h, AP.ki_h, AP.wn_h);
fprintf('  Va(pitch) kp = %9.4f  ki = %9.4f               wn = %6.3f\n', ...
        AP.kp_V2, AP.ki_V2, AP.wn_V2);
fprintf('  Va(throt) kp = %9.4f  ki = %9.4f               wn = %6.3f\n', ...
        AP.kp_V, AP.ki_V, AP.wn_V);

%% ======================================= open loop: elevator doublet
Ts = P.Ts;
T_end = 30;
N = round(T_end/Ts);
t = (0:N-1)'*Ts;

x = x_trim;
log_ol = zeros(N, 7);
for k = 1:N
    de = u_trim(1);
    if t(k) >= 2 && t(k) < 2.5
        de = de + deg2rad(2);
    elseif t(k) >= 2.5 && t(k) < 3
        de = de - deg2rad(2);
    end
    uk = [de; u_trim(2)];
    log_ol(k,:) = [x', uk'];
    x = rk4(x, uk, P, Ts);
end

%% ==================================== closed loop: +50 m altitude step
x = x_trim;
S = lon_autopilot_reset(u_trim(2));
T_end_cl = 120;
Ncl = round(T_end_cl/Ts);
tcl = (0:Ncl-1)'*Ts;

h_c  = x_trim(5) + 50;
Va_c = Va_star;

log_cl = zeros(Ncl, 8);
for k = 1:Ncl
    [uk, S, dbg] = lon_autopilot([h_c; Va_c], x, S, AP);
    log_cl(k,:) = [x', uk', dbg.theta_c];
    x = rk4(x, uk, P, Ts);
end

fprintf('\n--- Closed loop, %+g m altitude step ---\n', h_c - x_trim(5));
fprintf('  final h  = %8.3f m   (command %8.3f, error %+.3f)\n', ...
        log_cl(end,5), h_c, log_cl(end,5)-h_c);
fprintf('  final Va = %8.3f m/s (command %8.3f, error %+.3f)\n', ...
        hypot(log_cl(end,1),log_cl(end,2)), Va_c, ...
        hypot(log_cl(end,1),log_cl(end,2))-Va_c);

%% ============================================================= plots
figure('Name','Open loop: elevator doublet');
subplot(3,1,1); plot(t, hypot(log_ol(:,1),log_ol(:,2))); grid on;
ylabel('V_a [m/s]'); title('Open loop, 2 deg elevator doublet at t = 2 s');
subplot(3,1,2); plot(t, rad2deg(log_ol(:,4))); grid on; ylabel('\theta [deg]');
subplot(3,1,3); plot(t, log_ol(:,5)); grid on; ylabel('h [m]'); xlabel('t [s]');

figure('Name','Closed loop: altitude step');
subplot(4,1,1); plot(tcl, log_cl(:,5), tcl, h_c*ones(Ncl,1), '--'); grid on;
ylabel('h [m]'); legend('h','h_c','Location','southeast');
title(sprintf('Closed loop, %s', airframe));
subplot(4,1,2); plot(tcl, hypot(log_cl(:,1),log_cl(:,2)), tcl, Va_c*ones(Ncl,1), '--');
grid on; ylabel('V_a [m/s]');
subplot(4,1,3); plot(tcl, rad2deg(log_cl(:,4)), tcl, rad2deg(log_cl(:,8)), '--');
grid on; ylabel('\theta [deg]'); legend('\theta','\theta_c','Location','northeast');
subplot(4,1,4); plot(tcl, rad2deg(log_cl(:,6)), tcl, log_cl(:,7)*10); grid on;
ylabel('cmd'); legend('\delta_e [deg]','\delta_t \times10'); xlabel('t [s]');

fprintf('\nSaved lon_coef.mat. Run check_lon to validate against the book.\n');

%% ========================================================== helpers
function xn = rk4(x, u, P, dt)
k1 = f_lon(x,            u, P);
k2 = f_lon(x + dt/2*k1,  u, P);
k3 = f_lon(x + dt/2*k2,  u, P);
k4 = f_lon(x + dt*k3,    u, P);
xn = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

function show_tf(name, nd)
%SHOW_TF Print a {num,den} coefficient pair as a readable polynomial ratio.
% No Control System Toolbox needed. run_lon stores {num,den} in coef, so you can
% always do tf(coef.T_theta_delta_e{:}) yourself if you want an LTI object.
fprintf('  %s = %s / (%s)\n', name, poly2s(nd{1}), poly2s(nd{2}));
end

function s = poly2s(c)
%POLY2S Format a polynomial coefficient vector in descending powers of s.
c = c(:)';
n = numel(c) - 1;
s = '';
for k = 1:numel(c)
    p = n - k + 1;
    if c(k) == 0, continue; end
    if isempty(s)
        if c(k) == 1 && p > 0
            term = '';              % write s^2, not 1 s^2
        else
            term = sprintf('%.4g', c(k));
        end
    else
        term = sprintf(' %s %.4g', sgn(c(k)), abs(c(k)));
    end
    if p == 1
        term = [term ' s'];
    elseif p > 1
        term = sprintf('%s s^%d', term, p);
    end
    s = [s term];
end
s = strtrim(s);
if isempty(s), s = '0'; end
end

function c = sgn(v)
if v < 0, c = '-'; else, c = '+'; end
end
