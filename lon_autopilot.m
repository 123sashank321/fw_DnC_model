function [delta, S, dbg] = lon_autopilot(cmd, x, S, AP)
%LON_AUTOPILOT Longitudinal autopilot: four loops plus the altitude state machine.
%
%   [delta, S, dbg] = lon_autopilot(cmd, x, S, AP)
%
%   cmd = [h_c; Va_c]            commanded altitude m, airspeed m/s
%   x   = [u; w; q; theta; h]
%   S   = integrator state, from lon_autopilot_reset()
%
%   delta = [delta_e; delta_t]
%   dbg   = struct with theta_c and the active zone name
%
% Book section 6.4 and figure 6.14. Integrator state is carried in S rather than
% in persistent variables: check_lon runs the loop more than once, and persistent
% state would make the second run silently wrong.

h_c  = cmd(1);
Va_c = cmd(2);

q     = x(3);
theta = x(4);
h     = x(5);
Va    = hypot(x(1), x(2));

%% ----------------------------------------------- altitude state machine
if h <= AP.takeoff_h
    zone = 'takeoff';
elseif h < h_c - AP.altitude_zone
    zone = 'climb';
elseif h > h_c + AP.altitude_zone
    zone = 'descend';
else
    zone = 'hold';
end

switch zone
case 'takeoff'
    % Full throttle, fixed pitch. Airspeed is deliberately NOT regulated with
    % pitch here: pitching down just after take-off flies into the ground.
    delta_t = 1;
    theta_c = AP.takeoff_pitch;
    S.int_h  = 0;      % hold the outer integrators clear until the hold zone
    S.int_V2 = 0;

case 'climb'
    % Full throttle, airspeed held with pitch, which maximises climb rate and
    % avoids stall.
    delta_t = 1;
    [theta_c, S] = airspeed_from_pitch(Va_c, Va, S, AP);
    S.int_h = 0;

case 'descend'
    % Zero throttle, airspeed again held with pitch.
    delta_t = 0;
    [theta_c, S] = airspeed_from_pitch(Va_c, Va, S, AP);
    S.int_h = 0;

case 'hold'
    % Airspeed on throttle, altitude on pitch.
    [delta_t, S] = airspeed_from_throttle(Va_c, Va, S, AP);
    [theta_c, S] = altitude_from_pitch(h_c, h, S, AP);
    S.int_V2 = 0;
end

theta_c = sat(theta_c, -AP.theta_max, AP.theta_max);

%% ----------------------------------------- 6.4.1 pitch hold on elevator
% PD on pitch, with the rate term taken from q directly rather than by
% differentiating theta.
delta_e = AP.kp_theta*(theta_c - theta) - AP.kd_theta*q;
delta_e = sat(delta_e, -AP.delta_e_max, AP.delta_e_max);

delta = [delta_e; sat(delta_t, 0, 1)];

if nargout > 2
    dbg = struct('theta_c', theta_c, 'zone', zone, 'Va', Va);
end

end

%% ========================================================================
function [theta_c, S] = altitude_from_pitch(h_c, h, S, AP)
%ALTITUDE_FROM_PITCH PI on altitude error, output is a pitch command. Eq. 6.24-6.25.

e = h_c - h;
int_try = S.int_h + AP.Ts*e;
theta_c = AP.kp_h*e + AP.ki_h*int_try;

% Conditional integration: only accept the integrator update if the result is
% inside the pitch limit, so the integral cannot wind up against the saturation.
if abs(theta_c) <= AP.theta_max
    S.int_h = int_try;
else
    theta_c = sat(theta_c, -AP.theta_max, AP.theta_max);
end

end

%% ========================================================================
function [theta_c, S] = airspeed_from_pitch(Va_c, Va, S, AP)
%AIRSPEED_FROM_PITCH PI on airspeed error, output is a pitch command. Eq. 6.27-6.28.

e = Va_c - Va;
int_try = S.int_V2 + AP.Ts*e;
theta_c = AP.kp_V2*e + AP.ki_V2*int_try;

if abs(theta_c) <= AP.theta_max
    S.int_V2 = int_try;
else
    theta_c = sat(theta_c, -AP.theta_max, AP.theta_max);
end

end

%% ========================================================================
function [delta_t, S] = airspeed_from_throttle(Va_c, Va, S, AP)
%AIRSPEED_FROM_THROTTLE PI on airspeed error, output is throttle. Section 6.4.4.
%
% Trimmed about the cruise throttle: the PI supplies the increment, so the
% integrator starts from the trim value rather than from zero.

e = Va_c - Va;
int_try = S.int_V + AP.Ts*e;
delta_t = S.delta_t_trim + AP.kp_V*e + AP.ki_V*int_try;

if delta_t >= 0 && delta_t <= 1
    S.int_V = int_try;
else
    delta_t = sat(delta_t, 0, 1);
end

end

%% ========================================================================
function y = sat(u, lo, hi)
y = min(max(u, lo), hi);
end
