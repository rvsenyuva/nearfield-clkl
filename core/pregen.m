function [X_all,H_all,th_all,r_all,W_all,Y_all,R_all] = pregen(P, SNR_dB, use_exact)
%PREGEN  Pre-generate all Monte Carlo trial data before the estimation loop.
%
%  Pre-generating data into cell arrays has two benefits:
%    1. Enables "parfor" in the estimation loop without broadcast issues.
%    2. Keeps data generation (rand/randn calls) out of the timed section
%       in Fig. 6 so runtime measurements reflect estimator cost only.
%
%  [X_all,H_all,th_all,r_all,W_all,Y_all,R_all] = pregen(P, SNR_dB, use_exact)
%
%  All outputs are N_MC x 1 cell arrays.

N_MC = P.N_MC;
X_all = cell(N_MC,1);  H_all = cell(N_MC,1);
th_all = cell(N_MC,1); r_all  = cell(N_MC,1);
W_all  = cell(N_MC,1); Y_all  = cell(N_MC,1);
R_all  = cell(N_MC,1);

for mc = 1:N_MC
    [X_all{mc}, H_all{mc}, ~, th_all{mc}, r_all{mc}] = ...
        nf_gen_channel(P, SNR_dB, use_exact);
    [W_all{mc}, Y_all{mc}, R_all{mc}] = nf_hybrid_combiner(X_all{mc}, P);
end
end
