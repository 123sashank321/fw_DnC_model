function M = lon_modes(A_lon)
%LON_MODES Short-period and phugoid modes, exact and by the book's approximations.
%
%   M = lon_modes(A_lon)
%
% Book section 5.6. Returns a struct with, for each of the two modes, the
% eigenvalue pair, natural frequency and damping ratio, computed two ways:
%
%   .exact   eigenvalues of the leading 4x4 block of A_lon
%   .approx  the reduced-order formulas of section 5.6
%
% The 5th state is altitude, whose column in A_lon is identically zero, so it
% contributes only a pole at the origin and is excluded from the 4x4 block.
%
% Comparing the two is the point of section 5.6: the approximations are what the
% classical literature uses, and seeing where they diverge from the eigenvalues
% tells you when the decoupling assumptions break down.
%
% eig and roots only, so no Control System Toolbox is needed.

% Jacobian entries are exactly the stability derivatives of book table 5.2.
Xu = A_lon(1,1);  Xq = A_lon(1,3);
Zu = A_lon(2,1);  Zw = A_lon(2,2);  Zq = A_lon(2,3);
Mw = A_lon(3,2);  Mq = A_lon(3,3);

%% ------------------------------------------------------------- exact
lam = eig(A_lon(1:4,1:4));

% Sort by |imaginary part|: the short period is the fast, well-damped pair and
% the phugoid the slow, lightly damped one. Fall back to |lambda| if a pair has
% gone real, which happens for some airframes.
[~, idx] = sort(abs(imag(lam)), 'descend');
lam = lam(idx);
if abs(imag(lam(1))) < eps
    [~, idx] = sort(abs(lam), 'descend');
    lam = lam(idx);
end

M.exact.short   = pack(lam(1:2));
M.exact.phugoid = pack(lam(3:4));
M.exact.all     = lam;

%% ------------------------------------------- approximate, section 5.6
% Short period: u held constant, theta* = 0.
half = (Zw + Mq)/2;
disc = half^2 - Mq*Zw + Mw*Zq;
lam_sp = half + [1; -1]*sqrt(complex(disc));
M.approx.short = pack(lam_sp);

% Phugoid: alpha held constant, theta* = 0.
num  = Zu*Xq - Xu*Zq;
hp   = num/(2*Zq);
discp = hp^2 + 9.81*Zu/Zq;
lam_ph = -hp + [1; -1]*sqrt(complex(discp));
M.approx.phugoid = pack(lam_ph);

end

%% ========================================================================
function s = pack(lam)
%PACK Natural frequency and damping ratio from an eigenvalue pair.

lam = lam(:);
s.lambda = lam;
s.wn     = abs(lam(1));
if s.wn > 0
    s.zeta = -real(lam(1)) / s.wn;
else
    s.zeta = NaN;
end
s.period = NaN;
if abs(imag(lam(1))) > 0
    s.period = 2*pi/abs(imag(lam(1)));
end

end
