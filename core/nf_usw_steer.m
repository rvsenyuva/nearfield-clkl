function a = nf_usw_steer(theta, r, P)
%NF_USW_STEER  Exact uniform spherical-wave (USW) near-field steering vector.
%
%  a = nf_usw_steer(theta, r, P)
%
%  Implements eq. (3) of the blueprint (phase-only, amplitude uniform):
%
%      [a]_m = exp( -j*(2*pi/lambda) * ( ||p(theta,r) - s_m|| - r ) )
%
%  where the common-phase reference point is the array centre (r removed).
%
%  INPUTS
%    theta  : scalar angle of arrival [rad]
%    r      : scalar range from array centre [m]  (r > 0)
%    P      : parameter struct (nf_params)
%
%  OUTPUT
%    a      : M x 1 complex steering vector, ||a||^2 = M

M      = P.M;
lambda = P.lambda;
d_ant  = P.d_ant;

% Centred element indices  (eq. 1)
m_bar = ((0:M-1).' - (M-1)/2);   % M x 1

% Sensor 2-D positions: s_m = [m_bar*d_ant, 0]
sx = m_bar * d_ant;               % M x 1  (only x differs for a ULA)

% Path location in 2-D Cartesian (broadside = x-axis)
px = r * cos(theta);
py = r * sin(theta);

% Euclidean distance from path point to each sensor element
dist_m = sqrt((px - sx).^2 + py^2);   % M x 1

% Exact USW phase  (common-phase reference: distance r to array centre)
a = exp(-1j * (2*pi/lambda) * (dist_m - r));
end
