function a = nf_fresnel_steer(theta, u, P)
%NF_FRESNEL_STEER  Fresnel (second-order / chirp) near-field steering vector.
%
%  a = nf_fresnel_steer(theta, u, P)
%
%  Implements the chirp parameterisation of eq. (5):
%
%      [a]_m = exp( j*omega(theta)*m_bar  -  j*c(theta)*u*m_bar^2 )
%
%  with linear-phase coefficient   omega(theta) = 2*pi*d/lambda * cos(theta)
%  and curvature scaling           c(theta)     = pi*d^2/lambda * sin^2(theta)
%  so that kappa(theta,u) = c(theta)*u  and  u = 1/r.
%
%  INPUTS
%    theta : scalar angle [rad]
%    u     : scalar inverse-range  1/r  [1/m]
%    P     : parameter struct (nf_params)
%
%  OUTPUT
%    a     : M x 1 complex steering vector

M      = P.M;
lambda = P.lambda;
d_ant  = P.d_ant;

m_bar = ((0:M-1).' - (M-1)/2);   % M x 1

omega = (2*pi*d_ant/lambda) * cos(theta);          % linear phase slope
c_i   = (pi*d_ant^2/lambda) * sin(theta)^2;        % curvature scale  (eq. 10)
kappa = c_i * u;                                    % quadratic coefficient

a = exp(1j * omega * m_bar - 1j * kappa * m_bar.^2);
end


% Note: nf_fresnel_steer_du is in its own file (core/nf_fresnel_steer_du.m)
% so that MATLAB can resolve it from any calling context.
