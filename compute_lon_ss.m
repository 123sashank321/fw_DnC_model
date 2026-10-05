function [A_lon, B_lon] = compute_lon_ss(x, u, P, scheme, step)
%COMPUTE_LON_SS Linearise the longitudinal model by finite differences.
%
%   [A_lon, B_lon] = compute_lon_ss(x, u, P)
%   [A_lon, B_lon] = compute_lon_ss(x, u, P, 'forward', 0.01)
%
%   A_lon  5x5, rows/cols ordered [u, w, q, theta, h]
%   B_lon  5x2, cols ordered [delta_e, delta_t]
%
%   scheme  'central' (default) or 'forward'
%   step    perturbation. Default is a per-variable scaled step,
%           1e-6*max(1,|x_i|). A scalar here overrides it for every variable.
%
% Book Appendix F.4. Central differencing with a scaled step is used by default
% because it is O(h^2) accurate. The 'forward' / 0.01 combination exists only to
% reproduce mavsim's reference matrices, which were generated that way; that
% step is far too large for accuracy and is not a sensible default.
%
% Because the 5th state is h rather than p_d, the result comes out directly in
% the book's [u w q theta h] ordering with no row or column sign flip. mavsim
% integrates p_d and has to negate afterwards.

if nargin < 4 || isempty(scheme), scheme = 'central'; end
if nargin < 5, step = []; end

x = x(:);  u = u(:);
nx = numel(x);
nu = numel(u);

A_lon = zeros(nx, nx);
for i = 1:nx
    h = pert(step, x(i));
    A_lon(:,i) = diff_col(@(v) f_lon(v, u, P), x, i, h, scheme);
end

B_lon = zeros(nx, nu);
for j = 1:nu
    h = pert(step, u(j));
    B_lon(:,j) = diff_col(@(v) f_lon(x, v, P), u, j, h, scheme);
end

end

%% ========================================================================
function h = pert(step, value)
if isempty(step)
    h = 1e-6 * max(1, abs(value));
else
    h = step;
end
end

%% ========================================================================
function col = diff_col(fun, v, i, h, scheme)
%DIFF_COL One Jacobian column by perturbing element i of v.

switch lower(scheme)
case 'central'
    vp = v;  vp(i) = vp(i) + h;
    vm = v;  vm(i) = vm(i) - h;
    col = (fun(vp) - fun(vm)) / (2*h);
case 'forward'
    vp = v;  vp(i) = vp(i) + h;
    col = (fun(vp) - fun(v)) / h;
otherwise
    error('compute_lon_ss:scheme', ...
          'Unknown scheme "%s". Use ''central'' or ''forward''.', scheme);
end

end
