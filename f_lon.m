function [xdot, aux] = f_lon(x, u, P)
%F_LON Nonlinear longitudinal equations of motion, 5 states.
%
%   [xdot, aux] = f_lon(x, u, P)
%
%   x = [u; w; q; theta; h]     body velocities m/s, pitch rate rad/s,
%                               pitch angle rad, altitude m
%   u = [delta_e; delta_t]      elevator rad, throttle 0..1
%
%   aux returns Va, alpha, the coefficients and the component forces, for
%   plotting and for the transfer-function code that needs dT/dVa.
%
% This is the longitudinal reduction of book equations (5.3), (5.4), (5.6) and
% (5.11), with the lateral states forced to zero: v = p = r = phi = psi = 0.
% Under that reduction theta_dot = q exactly, and h is integrated directly
% instead of p_d, so no sign flip is needed anywhere downstream.
%
% Pure function, no persistent or global state: fzero, the finite-difference
% Jacobian and the RK4 loop all call it reentrantly.

ui    = x(1);
w     = x(2);
q     = x(3);
theta = x(4);

delta_e = u(1);
delta_t = u(2);

%% ------------------------------------------------------------- air data
Va    = sqrt(ui^2 + w^2);
alpha = atan2(w, ui);

% Floor the airspeed used in the damping denominators so the model stays
% finite at Va = 0. Same guard as ../NL_FW_Dynamics.
Va_safe = max(Va, 1e-3);
d_long  = P.c / (2*Va_safe);

qbar  = 0.5 * P.rho * Va^2;
qbarS = qbar * P.S;

%% ------------------------------------------------------ lift coefficient
CL_lin = P.C_L_0 + P.C_L_alpha*alpha;

if P.stall
    % Sigmoid blend from the linear model into a flat-plate model past alpha0
    % (book eq. 4.9-4.10). sigma -> 0 well inside the linear range.
    ex_m = exp(-P.M*(alpha - P.alpha0));
    ex_p = exp( P.M*(alpha + P.alpha0));
    sigma = (1 + ex_m + ex_p) / ((1 + ex_m)*(1 + ex_p));
    CL_flat = 2*sign(alpha)*sin(alpha)^2*cos(alpha);
    CL_static = (1 - sigma)*CL_lin + sigma*CL_flat;
else
    CL_static = CL_lin;
end

CL = CL_static + P.C_L_q*d_long*q + P.C_L_delta_e*delta_e;

%% ------------------------------------------------------ drag coefficient
% polyval covers both the book's linear-in-alpha form and the canard's
% quadratic polar; see lon_params.
CD = polyval(P.CD_poly, alpha) + P.C_D_q*d_long*q + P.C_D_delta_e*delta_e;

%% --------------------------------------------------------- pitch moment
% Note d_long uses the chord. Using the span here is a real bug present in the
% older canard Simulink block; see ../NL_FW_Dynamics/README.md.
Cm = P.C_m_0 + P.C_m_alpha*alpha + P.C_m_q*d_long*q + P.C_m_delta_e*delta_e;

%% ------------------------------------------------------------- forces
L  = qbarS * CL;
D  = qbarS * CD;
My = qbarS * P.c * Cm;

T = thrust(Va, delta_t, P);

% Rotate wind-frame lift and drag into body axes, add gravity and thrust.
ca = cos(alpha);  sa = sin(alpha);
fx = -P.mass*P.g*sin(theta) + (-D*ca + L*sa) + T;
fz =  P.mass*P.g*cos(theta) + (-D*sa - L*ca);

%% ------------------------------------------------------------ dynamics
xdot = [ -q*w + fx/P.mass
          q*ui + fz/P.mass
          My/P.Jy
          q
          ui*sin(theta) - w*cos(theta) ];

if nargout > 1
    aux = struct('Va',Va, 'alpha',alpha, 'CL',CL, 'CD',CD, 'Cm',Cm, ...
                 'L',L, 'D',D, 'My',My, 'T',T, 'fx',fx, 'fz',fz, 'qbar',qbar);
end

end

%% ========================================================================
function T = thrust(Va, delta_t, P)
%THRUST Propeller thrust along body x.

switch P.thrust_model

case 'prop_quad'
    % McLain addendum model. Match prop torque to motor torque to get the
    % propeller speed, then evaluate the thrust coefficient at that advance
    % ratio.
    v_in = P.V_max * delta_t;

    a = P.rho * P.D_prop^5 * P.C_Q0 / (4*pi^2);
    b = P.rho * P.D_prop^4 * P.C_Q1 * Va / (2*pi) + P.KQ^2 / P.R_motor;
    c = P.rho * P.D_prop^3 * P.C_Q2 * Va^2 - P.KQ*v_in/P.R_motor + P.KQ*P.i0;

    disc = b^2 - 4*a*c;
    if disc < 0
        T = 0;
        return
    end
    % -2c/(b+sqrt(disc)) rather than (-b+sqrt(disc))/(2a): the two agree
    % algebraically but this form does not lose precision when 4ac << b^2,
    % which is exactly the regime at cruise.
    Omega = -2*c / (b + sqrt(disc));
    if Omega <= 0
        T = 0;
        return
    end

    J   = 2*pi*Va / (Omega * P.D_prop);
    C_T = P.C_T2*J^2 + P.C_T1*J + P.C_T0;
    n   = Omega / (2*pi);
    T   = P.rho * n^2 * P.D_prop^4 * C_T;

case 'linear_lapse'
    T = P.T_static * (1 - Va/P.V_thrust_max) * delta_t;
    T = max(T, 0);

otherwise
    error('f_lon:thrust', 'Unknown thrust_model "%s".', P.thrust_model);
end

end
