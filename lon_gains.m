function AP = lon_gains(coef, P, tune)
%LON_GAINS Longitudinal autopilot gains by successive loop closure.
%
%   AP = lon_gains(coef, P)
%   AP = lon_gains(coef, P, tune)
%
% Book section 6.4. Four loops, designed inner to outer:
%
%   6.4.1  pitch attitude hold using elevator          -> kp_theta, kd_theta
%   6.4.2  altitude hold using commanded pitch         -> kp_h,  ki_h
%   6.4.3  airspeed hold using commanded pitch         -> kp_V2, ki_V2
%   6.4.4  airspeed hold using throttle                -> kp_V,  ki_V
%
% The outer loops are deliberately slower than the pitch loop by the bandwidth
% separation factors W_h and W_V2, which is what makes successive loop closure
% valid. Defaults put both at 10.
%
% tune fields (all optional): zeta_theta, W_h, zeta_h, W_V2, zeta_V2,
%                             wn_V, zeta_V, e_theta_max

if nargin < 3, tune = struct(); end

dflt = struct('zeta_theta', 0.707, ...
              'W_h',       10.0, 'zeta_h',   0.9, ...
              'W_V2',      10.0, 'zeta_V2',  0.9, ...
              'wn_V',       0.5, 'zeta_V',   0.9, ...
              'e_theta_max', deg2rad(10));
fn = fieldnames(dflt);
for k = 1:numel(fn)
    if ~isfield(tune, fn{k}) || isempty(tune.(fn{k}))
        tune.(fn{k}) = dflt.(fn{k});
    end
end

% Guard against designing gains for the wrong aircraft from a stale lon_coef.mat.
if isfield(coef, 'airframe') && ~strcmp(coef.airframe, P.airframe)
    error('lon_gains:airframeMismatch', ...
          ['coef was computed for "%s" but P is "%s". Re-run compute_lon_tf ' ...
           'or delete lon_coef.mat.'], coef.airframe, P.airframe);
end

Va = coef.Va_trim;
g  = P.g;

%% ------------------------------------------- 6.4.1 pitch hold (elevator)
% kp_theta is set so the elevator just saturates at the maximum anticipated
% pitch error. sign(a_theta3) is required for stability: a_theta3 follows
% C_m_delta_e, which is negative.
kp_theta = (P.delta_e_max / tune.e_theta_max) * sign(coef.a_theta3);
wn_theta = sqrt(coef.a_theta2 + kp_theta*coef.a_theta3);          % eq. 6.21
kd_theta = (2*tune.zeta_theta*wn_theta - coef.a_theta1) / coef.a_theta3;  % eq. 6.22
K_theta_DC = kp_theta*coef.a_theta3 / (coef.a_theta2 + kp_theta*coef.a_theta3);

%% ----------------------------------- 6.4.2 altitude hold (pitch command)
wn_h = wn_theta / tune.W_h;
ki_h = wn_h^2 / (K_theta_DC * Va);                               % eq. 6.24
kp_h = 2*tune.zeta_h*wn_h / (K_theta_DC * Va);                   % eq. 6.25

%% ----------------------------------- 6.4.3 airspeed hold (pitch command)
wn_V2 = wn_theta / tune.W_V2;
ki_V2 = -wn_V2^2 / (K_theta_DC * g);                             % eq. 6.27
kp_V2 = (coef.a_V1 - 2*tune.zeta_V2*wn_V2) / (K_theta_DC * g);   % eq. 6.28

%% --------------------------------------- 6.4.4 airspeed hold (throttle)
wn_V = tune.wn_V;
ki_V = wn_V^2 / coef.a_V2;
kp_V = (2*tune.zeta_V*wn_V - coef.a_V1) / coef.a_V2;

%% ------------------------------------------------------------- assemble
AP = struct();
AP.airframe = P.airframe;
AP.Ts       = P.Ts;
AP.Va0      = Va;
AP.g        = g;

AP.kp_theta   = kp_theta;
AP.kd_theta   = kd_theta;
AP.K_theta_DC = K_theta_DC;
AP.wn_theta   = wn_theta;

AP.kp_h = kp_h;  AP.ki_h = ki_h;  AP.wn_h  = wn_h;
AP.kp_V2 = kp_V2; AP.ki_V2 = ki_V2; AP.wn_V2 = wn_V2;
AP.kp_V = kp_V;  AP.ki_V = ki_V;  AP.wn_V  = wn_V;

% Limits and the altitude state machine thresholds
AP.delta_e_max   = P.delta_e_max;
AP.theta_max     = P.theta_max;
AP.altitude_zone = 10;             % m, half-width of the altitude-hold band
AP.takeoff_h     = 5;              % m, end of the take-off zone
AP.takeoff_pitch = deg2rad(10);

if wn_h >= wn_theta || wn_V2 >= wn_theta
    warning('lon_gains:bandwidth', ...
            ['Outer-loop bandwidth is not below the pitch loop ' ...
             '(wn_theta=%.3f, wn_h=%.3f, wn_V2=%.3f). Successive loop ' ...
             'closure assumes separation; increase W_h / W_V2.'], ...
            wn_theta, wn_h, wn_V2);
end

end
