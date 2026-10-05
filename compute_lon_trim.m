function [x_trim, u_trim, info] = compute_lon_trim(Va, gamma, P)
%COMPUTE_LON_TRIM Wings-level longitudinal trim at airspeed Va, climb angle gamma.
%
%   [x_trim, u_trim, info] = compute_lon_trim(Va, gamma, P)
%
%   x_trim = [u; w; q; theta; h]      u_trim = [delta_e; delta_t]
%   info   = struct with alpha, theta, residual, and the solver exit flags
%
% Book Appendix F.2, reduced to the longitudinal plane. The problem is fully
% TRIANGULAR, so no multivariable optimiser is needed:
%
%   1. alpha fixes the kinematics:  theta = alpha + gamma, q = 0,
%      u = Va*cos(alpha), w = Va*sin(alpha)
%   2. delta_e is then CLOSED FORM from q_dot = 0 with q = 0:
%         delta_e = -(C_m_0 + C_m_alpha*alpha) / C_m_delta_e
%   3. alpha is a scalar root-find on w_dot = 0. Thrust acts along body x and
%      so does not appear in f_z, which is what makes this step independent of
%      delta_t.
%   4. delta_t is a scalar root-find on u_dot = 0, done last.
%
% h_dot needs no equation: u*sin(theta) - w*cos(theta) = Va*sin(gamma) holds
% identically once theta = alpha + gamma, so the climb rate is automatic.
%
% Uses fzero only, so there is no Optimization Toolbox dependency. mavsim uses
% fmincon on a 17-element vector because it trims all 12 states at once.

if nargin < 2 || isempty(gamma), gamma = 0; end

%% ------------------------------- step 3: alpha from the w_dot = 0 residual
bracket = [-0.2, 0.3];          % rad; physically meaningful, ~-11 to +17 deg
res = @(a) wdot_residual(a, Va, gamma, P);

ra = res(bracket(1));
rb = res(bracket(2));
if ~isfinite(ra) || ~isfinite(rb)
    error('compute_lon_trim:residual', ...
          'w_dot residual is not finite at the bracket ends.');
end
if sign(ra) == sign(rb)
    error('compute_lon_trim:noTrim', ...
          ['No longitudinal trim for Va = %.4g m/s, gamma = %.4g rad: the ' ...
           'w_dot residual does not change sign on alpha in [%.3g, %.3g] ' ...
           '(residuals %.4g and %.4g). The aircraft cannot hold this ' ...
           'condition.'], Va, gamma, bracket(1), bracket(2), ra, rb);
end

[alpha, ~, flag_a] = fzero(res, bracket);
if flag_a ~= 1
    error('compute_lon_trim:alphaFzero', 'fzero on alpha failed, flag %d.', flag_a);
end

%% --------------------------------------- steps 1-2: kinematics and elevator
[x_trim, delta_e] = state_from_alpha(alpha, Va, gamma, P);

%% ----------------------------- step 4: delta_t from the u_dot = 0 residual
% Not clamped to [0,1] inside f_lon: for 'prop_quad' the thrust coefficient
% polynomial extrapolates to negative thrust at low throttle / high advance
% ratio, which is unphysical but keeps the residual smooth and sign-changing
% so fzero stays well posed. The returned delta_t is checked below.
ures = @(dt) subsref_1(f_lon(x_trim, [delta_e; dt], P));

u0 = ures(0);
u1 = ures(1);
if sign(u0) == sign(u1)
    error('compute_lon_trim:noThrottle', ...
          ['No throttle in [0,1] trims u_dot at Va = %.4g m/s, gamma = ' ...
           '%.4g rad (u_dot = %.4g at delta_t = 0 and %.4g at 1). The ' ...
           'propulsion model cannot sustain this airspeed.'], ...
          Va, gamma, u0, u1);
end

[delta_t, ~, flag_t] = fzero(ures, [0, 1]);
if flag_t ~= 1
    error('compute_lon_trim:throttleFzero', 'fzero on delta_t failed, flag %d.', flag_t);
end

u_trim = [delta_e; delta_t];

%% -------------------------------------------------------------- report out
xdot = f_lon(x_trim, u_trim, P);
info = struct('alpha', alpha, ...
              'theta', x_trim(4), ...
              'Va', Va, ...
              'gamma', gamma, ...
              'residual', norm(xdot(1:4)), ...
              'hdot', xdot(5), ...
              'hdot_expected', Va*sin(gamma));

end

%% ========================================================================
function [x, delta_e] = state_from_alpha(alpha, Va, gamma, P)
%STATE_FROM_ALPHA Kinematics and the closed-form elevator for a given alpha.

theta = alpha + gamma;
x = [ Va*cos(alpha)
      Va*sin(alpha)
      0
      theta
      P.h0 ];

% q_dot = 0 with q = 0  =>  C_m_0 + C_m_alpha*alpha + C_m_delta_e*delta_e = 0
delta_e = -(P.C_m_0 + P.C_m_alpha*alpha) / P.C_m_delta_e;

end

%% ========================================================================
function r = wdot_residual(alpha, Va, gamma, P)
%WDOT_RESIDUAL w_dot at the alpha-implied state. Independent of throttle.

[x, delta_e] = state_from_alpha(alpha, Va, gamma, P);
xdot = f_lon(x, [delta_e; 0], P);
r = xdot(2);

end

%% ========================================================================
function v = subsref_1(xdot)
%SUBSREF_1 First element of a vector. MATLAB forbids f(...)(1) inline.
v = xdot(1);
end
