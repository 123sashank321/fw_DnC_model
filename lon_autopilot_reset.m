function S = lon_autopilot_reset(delta_t_trim)
%LON_AUTOPILOT_RESET Fresh integrator state for lon_autopilot.
%
%   S = lon_autopilot_reset(delta_t_trim)
%
% delta_t_trim is the cruise throttle the airspeed-on-throttle loop works around.
% Defaults to 0.5 if omitted.

if nargin < 1 || isempty(delta_t_trim)
    delta_t_trim = 0.5;
end

S = struct('int_h', 0, 'int_V2', 0, 'int_V', 0, 'delta_t_trim', delta_t_trim);

end
