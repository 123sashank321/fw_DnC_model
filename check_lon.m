%% check_lon.m -- the self-check for the longitudinal chain
%
% Four groups, in increasing scope. Groups 1, 2 and 4 are reference-free and run
% for BOTH airframes. Group 3 compares against the published mavsim numbers and
% so runs for the Aerosonde only.
%
% Reference values are from
%   mavsim_public/mavsim_python/launch_files/chap05/chap5_check.py
% at Va = 25, gamma = 0, generated with the PYTHON parameter set.
%
% IMPORTANT, and the reason group 3 is split in two: the published trim point is
% NOT a force equilibrium for the parameters that ship with it. f_lon there gives
% u_dot = -0.85 m/s^2, because the propeller model makes only about 0.95 N of
% thrust at delta_t = 0.6768 while the drag is 10.3 N. So we do NOT compare our
% trim against theirs. Instead we linearise at THEIR point to compare A_lon and
% B_lon apples to apples, and separately require OUR trim to be a true
% equilibrium. Both checks are then tight.
%
% Run:  check_lon

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);

ok = true(0);
names = {};

fprintf('\n================ check_lon ================\n');

%% ============================ group 1: trim is a true equilibrium
fprintf('\n[1] Trim residual, both airframes, level and climbing\n');
for af = {'aerosonde','canard'}
    P = lon_params(af{1});
    for gam = [0, 0.1]
        [~, ut, info] = compute_lon_trim(P.Va0, gam, P);
        r = info.residual;
        % Climb rate must come out right on its own; that is what catches a
        % gravity or theta sign error which cancels at gamma = 0.
        dh = abs(info.hdot - P.Va0*sin(gam));
        [ok, names] = rec(ok, names, ...
            sprintf('%s trim resid (gamma=%.2f)', af{1}, gam), r, 1e-9);
        [ok, names] = rec(ok, names, ...
            sprintf('%s climb rate  (gamma=%.2f)', af{1}, gam), dh, 1e-9);
        fprintf('    %-10s gamma=%.2f  |f|=%.2e  alpha=%7.4f  de=%8.5f  dt=%6.4f\n', ...
                af{1}, gam, r, info.alpha, ut(1), ut(2));
    end
end

%% ======== group 2: Jacobian agrees with the analytic a-coefficients
fprintf('\n[2] A_lon / B_lon against the analytic coefficients, both airframes\n');
for af = {'aerosonde','canard'}
    P = lon_params(af{1});
    [xt, ut] = compute_lon_trim(P.Va0, 0, P);
    [A, B]   = compute_lon_ss(xt, ut, P);
    c        = compute_lon_tf(xt, ut, P);

    % Two independent derivations of the q_dot row: numerical Jacobian versus
    % the closed-form section 5.4 coefficients. Agreement is a far stronger
    % statement than either matching a stored number.
    [ok, names] = rec(ok, names, sprintf('%s A(3,3) = -a_theta1', af{1}), ...
                      relerr(A(3,3), -c.a_theta1), 1e-6);
    [ok, names] = rec(ok, names, sprintf('%s B(3,1) =  a_theta3', af{1}), ...
                      relerr(B(3,1),  c.a_theta3), 1e-6);
    % a_theta2 is defined against theta, the Jacobian against w; the two relate
    % only through alpha ~ w/Va, so this one is looser by construction.
    [ok, names] = rec(ok, names, sprintf('%s A(3,2) = -a_theta2/Va', af{1}), ...
                      relerr(A(3,2), -c.a_theta2/c.Va_trim), 1e-2);
    % theta_dot = q exactly.
    [ok, names] = rec(ok, names, sprintf('%s A(4,3) = 1', af{1}), ...
                      abs(A(4,3)-1), 1e-9);
    % Altitude feeds back into nothing, so its column is identically zero.
    [ok, names] = rec(ok, names, sprintf('%s A(:,5) = 0', af{1}), ...
                      norm(A(:,5)), 1e-9);
    fprintf('    %-10s A(3,3)=%11.6f  -a_theta1=%11.6f\n', af{1}, A(3,3), -c.a_theta1);
end

%% ================= group 3: against the published mavsim numbers
fprintf('\n[3] Against chap5_check.py, Aerosonde only\n');
P = lon_params('aerosonde');

% --- 3a: coefficients that depend only on Va and the airframe
[xt, ut] = compute_lon_trim(25, 0, P);
c = compute_lon_tf(xt, ut, P);
[ok, names] = rec(ok, names, 'a_theta1', relerr(c.a_theta1,   5.294738359662443), 1e-5);
[ok, names] = rec(ok, names, 'a_theta2', relerr(c.a_theta2,  99.94742395724161 ), 1e-4);
[ok, names] = rec(ok, names, 'a_theta3', relerr(c.a_theta3, -36.11239040790846 ), 1e-5);
[ok, names] = rec(ok, names, 'a_V3',     relerr(c.a_V3,       9.809999999999892), 1e-6);
fprintf('    a_theta1 %12.6f  vs %12.6f\n', c.a_theta1,  5.294738359662443);
fprintf('    a_theta3 %12.6f  vs %12.6f\n', c.a_theta3, -36.11239040790846);

% --- 3b: quantities that depend on the trim point, evaluated at THEIR point
x_ref = [24.968743; 1.249755; 0; 0.05001119284259148; 100];
u_ref = [-0.124778; 0.676752];
cR = compute_lon_tf(x_ref, u_ref, P);
[ok, names] = rec(ok, names, 'a_V1   @ ref pt', relerr(cR.a_V1,    0.2888454121899283), 1e-3);
[ok, names] = rec(ok, names, 'dT_dVa @ ref pt', relerr(cR.dT_dVa, -2.3521982958853327), 1e-3);
% a_V2 / dT_ddelta_t are looser: thrust is strongly curved in throttle, so the
% reference's own differencing step shows up here. Ours is the more accurate value.
[ok, names] = rec(ok, names, 'a_V2   @ ref pt', relerr(cR.a_V2,    8.20722086381344 ), 2e-2);

% --- 3c: the matrices, differenced the way the reference was
A_ref = [-0.27486079, 0.49868179, -1.21983882, -9.79511927, 0;
         -0.56234374,-4.49810469, 24.37105023, -0.53938541, 0;
          0.19993539,-3.99297865, -5.29473836,  0         , 0;
          0         , 0         ,  0.99997406,  0         , 0;
          0.04999035,-0.99874970,  0         , 24.99958361, 0];
B_ref = [ -0.13840016,  8.20722086;
          -2.58618345,  0;
         -36.11239041,  0;
           0         ,  0;
           0         ,  0];
[Af, Bf] = compute_lon_ss(x_ref, u_ref, P, 'forward', 0.01);
dA = max(abs(Af(:)-A_ref(:))) / max(abs(A_ref(:)));
dB = max(abs(Bf(:)-B_ref(:))) / max(abs(B_ref(:)));
[ok, names] = rec(ok, names, 'A_lon @ ref pt (fwd/0.01)', dA, 1e-3);
[ok, names] = rec(ok, names, 'B_lon @ ref pt (fwd/0.01)', dB, 1e-3);
fprintf('    A_lon scaled max dev %.3e,  B_lon %.3e\n', dA, dB);

% Record, do not assert: the published trim is not an equilibrium.
xd = f_lon(x_ref, u_ref, P);
fprintf('    NOTE published trim is not an equilibrium: |f| = %.4g (u_dot = %.4g)\n', ...
        norm(xd(1:4)), xd(1));

%% ===================== group 4: closed loop, chapter 6, both airframes
fprintf('\n[4] Closed-loop altitude step, both airframes\n');
for af = {'aerosonde','canard'}
    P = lon_params(af{1});
    [xt, ut] = compute_lon_trim(P.Va0, 0, P);
    c  = compute_lon_tf(xt, ut, P);
    AP = lon_gains(c, P);

    x = xt;  S = lon_autopilot_reset(ut(2));
    h_c = xt(5) + 10;  Va_c = P.Va0;
    N = round(120/P.Ts);
    theta_pk = 0;  theta_c_pk = 0;  dt_bad = false;
    for k = 1:N
        [uk, S, dbg] = lon_autopilot([h_c; Va_c], x, S, AP);
        theta_pk   = max(theta_pk,   abs(x(4)));
        theta_c_pk = max(theta_c_pk, abs(dbg.theta_c));
        if uk(2) < -1e-12 || uk(2) > 1+1e-12, dt_bad = true; end
        x = rk4c(x, uk, P, P.Ts);
    end
    h_err  = abs(x(5) - h_c);
    Va_err = abs(hypot(x(1),x(2)) - Va_c);

    [ok, names] = rec(ok, names, sprintf('%s h -> h_c',     af{1}), h_err,  0.5);
    [ok, names] = rec(ok, names, sprintf('%s Va held',      af{1}), Va_err, 2.0);
    % The autopilot limits the pitch COMMAND, so that is what gets the tight
    % check. The pitch state is a damped second-order response to that command
    % and so legitimately overshoots it; allow 25% and use it only as a
    % divergence guard.
    [ok, names] = rec(ok, names, sprintf('%s theta_c limit', af{1}), ...
                      max(0, theta_c_pk - AP.theta_max), 1e-9);
    [ok, names] = rec(ok, names, sprintf('%s theta bounded', af{1}), ...
                      max(0, theta_pk - 1.25*AP.theta_max), 1e-9);
    [ok, names] = rec(ok, names, sprintf('%s throttle 0..1',af{1}), double(dt_bad), 0.5);
    fprintf('    %-10s h err %6.3f m  Va err %6.3f m/s  peak |theta| %5.2f deg (cmd %5.2f)\n', ...
            af{1}, h_err, Va_err, rad2deg(theta_pk), rad2deg(theta_c_pk));
end

%% ======================================================== verdict
fprintf('\n---------------- results ----------------\n');
for i = 1:numel(names)
    fprintf('  %-34s %s\n', names{i}, tern(ok(i),'PASS','FAIL'));
end
nf = sum(~ok);
fprintf('-----------------------------------------\n');
if nf == 0
    fprintf('  ALL %d CHECKS PASSED\n\n', numel(ok));
else
    fprintf('  %d of %d FAILED\n\n', nf, numel(ok));
end
assert(nf == 0, 'check_lon: %d of %d checks failed.', nf, numel(ok));

%% ======================================================== helpers
function [ok, names] = rec(ok, names, name, value, tol)
ok(end+1)    = isfinite(value) && value <= tol;
names{end+1} = name;
end

function e = relerr(a, b)
if abs(b) > 1e-12
    e = abs(a-b)/abs(b);
else
    e = abs(a-b);
end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function xn = rk4c(x, u, P, dt)
k1 = f_lon(x,            u, P);
k2 = f_lon(x + dt/2*k1,  u, P);
k3 = f_lon(x + dt/2*k2,  u, P);
k4 = f_lon(x + dt*k3,    u, P);
xn = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end
