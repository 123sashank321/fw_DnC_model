function P = lon_params(airframe)
%LON_PARAMS Longitudinal parameters for one airframe, in canonical book names.
%
%   P = lon_params('aerosonde')   Beard & McLain Aerosonde
%   P = lon_params('canard')      the canard config used by ../NL_FW_Dynamics
%
% Both branches emit the SAME flat struct with the book's field names, SI units,
% angles in radians. Only two model tags differ between airframes: P.stall and
% P.thrust_model. Drag needs no tag because P.CD_poly is a polyval coefficient
% vector, and the book's "linear in alpha" drag is just a degenerate quadratic.
%
% Reference: Beard & McLain, "Small Unmanned Aircraft: Theory and Practice",
%            ch. 3-4 for the model, Appendix E for the Aerosonde data.

if nargin < 1 || isempty(airframe)
    airframe = 'aerosonde';
end

switch lower(airframe)

%% ------------------------------------------------------------- Aerosonde
case 'aerosonde'
    % !! THESE ARE THE PYTHON VALUES, NOT THE MATLAB ONES. !!
    % mavsim ships two disagreeing copies of the Aerosonde data. The chapter-5
    % reference numbers in mavsim_python/launch_files/chap05/chap5_check.py were
    % generated with the Python set. Verified: C_D_0 = 0.0424 and C_D_alpha =
    % 0.132 reproduce a_V1 = 0.2888454 to 7 digits; the MATLAB set (0.043 and
    % 0.030) does not. Do not "correct" these back to the MATLAB file's values.
    P.airframe = 'aerosonde';

    P.mass = 11.0;
    P.Jy   = 1.135;
    P.g    = 9.81;
    P.rho  = 1.2682;
    P.S    = 0.55;
    P.c    = 0.18994;        % chord (Python value; MATLAB file rounds to 0.19)
    P.b    = 2.8956;         % span, carried for reference only

    P.C_L_0       =  0.23;
    P.C_L_alpha   =  5.61;
    P.C_L_q       =  7.95;
    P.C_L_delta_e =  0.13;

    % CD = polyval(CD_poly, alpha) = 0.132*alpha + 0.0424
    P.CD_poly     = [0, 0.132, 0.0424];
    P.C_D_q       =  0.0;
    P.C_D_delta_e =  0.0135;

    P.C_m_0       =  0.0135;
    P.C_m_alpha   = -2.74;
    P.C_m_q       = -38.21;
    P.C_m_delta_e = -0.99;

    % sigmoid blend into a flat-plate post-stall model (book eq. 4.9-4.10)
    P.stall  = true;
    P.M      = 50.0;
    P.alpha0 = 0.47;

    % Propeller / motor, McLain addendum model. The legacy
    % (S_prop, k_motor, C_prop) model is NOT implemented: it cannot reproduce
    % the reference dT_dVa (it gives -6.43 against the required -2.35).
    P.thrust_model = 'prop_quad';
    P.D_prop  = 20*0.0254;                  % 0.508 m
    P.KV      = (1/145)*60/(2*pi);          % back-emf constant, V-s/rad
    P.KQ      = P.KV;                       % torque constant, N-m/A
    P.R_motor = 0.042;                      % ohm
    P.i0      = 1.5;                        % no-load current, A
    P.V_max   = 3.7*12;                     % 12 cells, 44.4 V
    P.C_Q2 = -0.01664;  P.C_Q1 = 0.004970;  P.C_Q0 = 0.005230;
    P.C_T2 = -0.1079;   P.C_T1 = -0.06044;  P.C_T0 = 0.09357;

    P.Va0 = 25;
    P.h0  = 100;

%% ---------------------------------------------------------------- Canard
case 'canard'
    % Values transcribed from ../NL_FW_Dynamics/init_params_fw.m, which in turn
    % took them from the canard configuration in
    % ../Dynamic Model/UAV_Non_linear_Dynamics/init_params.m (the linear
    % stability-derivative set). Inlined rather than loaded so this folder is
    % self-contained; if init_params_fw.m changes, update here too.
    P.airframe = 'canard';

    P.mass = 2.4;
    P.Jy   = 0.08381;
    P.g    = 9.814;
    P.rho  = 1.225;
    P.S    = 0.190;          % aero.S_w
    P.c    = 0.199;          % aero.c_bar
    P.b    = 1.0;

    P.C_L_0       = 0.3494;  % aero.CLo
    P.C_L_alpha   = 4.5106;  % aero.CL_a
    P.C_L_q       = 5.1889;
    P.C_L_delta_e = 0.43;

    % CD = 2.8252*alpha^2 - 0.0137*alpha + 0.0461  (quadratic polar in alpha)
    P.CD_poly     = [2.8252, -0.0137, 0.0461];
    P.C_D_q       = 0.0;
    P.C_D_delta_e = 0.06;

    P.C_m_0       =  0.0121;
    P.C_m_alpha   = -0.3263;
    P.C_m_q       = -5.451;
    P.C_m_delta_e = -1.28;

    % No stall model: alpha0 / M / epsilon were never identified for this
    % airframe. Do not invent them, and do not try to encode "no stall" as
    % alpha0 = inf -- the sigmoid returns NaN.
    P.stall = false;

    % Velocity-lapsed thrust, clamped non-negative.
    % P.T_lever_arm from init_params_fw.m is deliberately NOT applied: thrust
    % produces no pitching moment here, matching the book's model.
    P.thrust_model  = 'linear_lapse';
    P.T_static      = 1.767*9.814*2;        % both engines at Va = 0, 34.68 N
    P.V_thrust_max  = 70;

    P.Va0 = 25;
    P.h0  = 100;

otherwise
    error('lon_params:airframe', ...
          'Unknown airframe "%s". Use ''aerosonde'' or ''canard''.', airframe);
end

%% ------------------------------------------- shared simulation / AP limits
P.Ts          = 0.01;              % s
P.delta_e_max = deg2rad(45);
P.theta_max   = deg2rad(30);

end
