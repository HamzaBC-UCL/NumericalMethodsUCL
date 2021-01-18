function [VcMC,VpMC,cputime_MC] = mcBSM(nblocks,muRN,r,sigma,nsample,T,S0,K)

% Function that produces the monte carlo prices of a call and put option
% + CPU time

tic;
VcMCb = zeros(nblocks,1);
VpMCb = zeros(nblocks,1);
for i = 1:nblocks

    % Arithmetic Brownian motion X = log S at time T
    X = muRN*T + sigma*randn(1,nsample)*sqrt(T);

    % Geometric Brownian motion at time T
    S = S0*exp(X);

    % Discounted expected payoff
    VcMCb(i) = exp(-r*T)*mean(max(S-K,0));
    VpMCb(i) = exp(-r*T)*mean(max(K-S,0));

end
VcMC = mean(VcMCb);
VpMC = mean(VpMCb);
scMC = sqrt(var(VcMCb)/nblocks);
spMC = sqrt(var(VpMCb)/nblocks);
cputime_MC = toc;
fprintf('%20s%14.10f%14.10f%14.10f\n','Monte Carlo',VcMC,VpMC,cputime_MC)
fprintf('%20s%14.10f%14.10f\n','Monte Carlo stdev',scMC,spMC)

end



