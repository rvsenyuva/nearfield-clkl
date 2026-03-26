function [da_du, c_i] = nf_fresnel_steer_du(theta, a, P)
%NF_FRESNEL_STEER_DU  Derivative of Fresnel steering vector w.r.t. u = 1/r.
%
%  [da_du, c_i] = nf_fresnel_steer_du(theta, a, P)
%
%  From blueprint eq. (19):
%      d a_i / d u_i  =  -j * c_i * ( m_bar.^2  .*  a_i )
%
%  where  c_i = (pi * d_ant^2 / lambda) * sin^2(theta)
%  is the curvature scaling constant (blueprint eq. 10).
%
%  INPUTS
%    theta  : scalar angle [rad]
%    a      : M x 1  current Fresnel steering vector  a_i(u_i)
%             (computed beforehand via nf_fresnel_steer)
%    P      : parameter struct (nf_params)
%
%  OUTPUTS
%    da_du  : M x 1  complex derivative vector
%    c_i    : scalar curvature scale constant (reusable in gradient eqs.)

M      = P.M;
lambda = P.lambda;
d_ant  = P.d_ant;

m_bar = ((0:M-1).' - (M-1)/2);               % M x 1  centred indices
c_i   = (pi * d_ant^2 / lambda) * sin(theta)^2;

da_du = -1j * c_i * (m_bar.^2 .* a);
end
