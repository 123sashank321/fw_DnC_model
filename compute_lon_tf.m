function coef = compute_lon_tf(x, u, P)
%COMPUTE_LON_TF Longitudinal transfer functions and their a-coefficients.
%
%   coef = compute_lon_tf(x_trim, u_trim, P)
%
% Book section 5.4, longitudinal part. Returns a struct holding
%
%   a_theta1 a_theta2 a_theta3      pitch from elevator      (eq. 5.29)
%   a_V1 a_V2 a_V3                  airspeed                 (eq. 5.35-5.36)
%   dT_dVa dT_ddelta_t              thrust sensitivities
%   Va_trim alpha_trim theta_trim delta_e_trim delta_t_trim
%
% and the five transfer functions as {num, den} coefficient-vector pairs:
%
%   T_theta_delta_e = a_theta3 / (s^2 + a_theta1*s + a_theta2)
%   T_h_theta       = Va / s
%   T_h_Va          = theta / s
%   T_Va_delta_t    = a_V2 / (s + a_V1)
%   T_Va_theta      = -a_V3 / (s + a_V1)
%
% Coefficient vectors, not tf objects: Control System Toolbox is never required
% here. run_lon builds tf objects for display only if the toolbox is present.

x = x(:);  u = u(:);
[~, aux] = f_lon(x, u, P);

Va      = aux.Va;
alpha   = aux.alpha;
theta   = x(4);
delta_e = u(1);
delta_t = u(2);
gamma   = theta - alpha;

%% ------------------------------------------------ pitch from elevator
% Common factor rho*Va^2*c*S/(2*Jy).
k_m = P.rho * Va^2 * P.c * P.S / (2*P.Jy);

a_theta1 = -k_m * P.C_m_q * P.c / (2*Va);
a_theta2 = -k_m * P.C_m_alpha;
a_theta3 =  k_m * P.C_m_delta_e;

%% ----------------------------------------------------------- airspeed
% Static drag coefficient at trim: the alpha polynomial plus the elevator term.
% The q term drops out because q = 0 in trim.
CD_trim = polyval(P.CD_poly, alpha) + P.C_D_delta_e*delta_e;

dT_dVa      = dthrust(x, u, P, 'Va');
dT_ddelta_t = dthrust(x, u, P, 'delta_t');

% Book eq. 5.35 writes the propulsion term for the legacy prop model as
% +rho*S_prop*C_prop*Va/m, which is exactly -dT_dVa/m. Written that way here so
% it holds for either thrust model.
a_V1 = P.rho * Va * P.S / P.mass * CD_trim - dT_dVa / P.mass;
a_V2 = dT_ddelta_t / P.mass;
a_V3 = P.g * cos(gamma);

%% ----------------------------------------------------------- assemble
coef = struct();
coef.airframe = P.airframe;

coef.Va_trim      = Va;
coef.alpha_trim   = alpha;
coef.theta_trim   = theta;
coef.gamma_trim   = gamma;
coef.delta_e_trim = delta_e;
coef.delta_t_trim = delta_t;

coef.a_theta1 = a_theta1;
coef.a_theta2 = a_theta2;
coef.a_theta3 = a_theta3;
coef.a_V1     = a_V1;
coef.a_V2     = a_V2;
coef.a_V3     = a_V3;

coef.dT_dVa      = dT_dVa;
coef.dT_ddelta_t = dT_ddelta_t;

coef.T_theta_delta_e = {a_theta3,  [1, a_theta1, a_theta2]};
coef.T_h_theta       = {Va,        [1, 0]};
coef.T_h_Va          = {theta,     [1, 0]};
coef.T_Va_delta_t    = {a_V2,      [1, a_V1]};
coef.T_Va_theta      = {-a_V3,     [1, a_V1]};

end

%% ========================================================================
function d = dthrust(x, u, P, wrt)
%DTHRUST Central-difference sensitivity of thrust to Va or to delta_t.
%
% Differentiating f_lon's reported thrust rather than re-deriving it keeps the
% two thrust models in one place. Perturbing Va is done by scaling u and w so
% alpha is held fixed, which is what the partial derivative means here.

switch wrt
case 'Va'
    Va = hypot(x(1), x(2));
    h  = 1e-6 * max(1, Va);
    d  = (T_at_Va(x, u, P, Va + h) - T_at_Va(x, u, P, Va - h)) / (2*h);
case 'delta_t'
    h  = 1e-6 * max(1, abs(u(2)));
    up = u;  up(2) = up(2) + h;
    um = u;  um(2) = um(2) - h;
    [~, ap] = f_lon(x, up, P);
    [~, am] = f_lon(x, um, P);
    d = (ap.T - am.T) / (2*h);
otherwise
    error('compute_lon_tf:dthrust', 'wrt must be ''Va'' or ''delta_t''.');
end

end

%% ========================================================================
function T = T_at_Va(x, u, P, Va_new)
%T_AT_VA Thrust with the velocity vector rescaled to Va_new, alpha unchanged.

Va = hypot(x(1), x(2));
xs = x;
xs(1) = x(1) * Va_new / Va;
xs(2) = x(2) * Va_new / Va;
[~, aux] = f_lon(xs, u, P);
T = aux.T;

end
