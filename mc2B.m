function [VcMC2B,VpMC2B,cputime_MC2B] = mc2B(nblocks,muRN2A,r,sigma,nsample,T,S0,K,lamda,muj,sigmaj)

% Function that produces merton jump monte carlo prices of a call and put option
% + CPU time

tic;
VcMCb = zeros(nblocks,1);
VpMCb = zeros(nblocks,1);
for i = 1:nblocks

    % Arithmetic Brownian motion 
    ADx = muRN2A*T + sigma*randn(1,nsample)*sqrt(T);
    
    % increments of NCPP
    dN = poissrnd(lamda*T,[1,nsample]);
    
    JDx = muj*dN + sigmaj*sqrt(dN).*randn(1,nsample); 
    
    % sum of increments NCPP & ABM
    X = JDx + ADx;
    
    % Geometric Brownian motion at time T
    S = S0*exp(X);

    % Discounted expected payoff
    VcMCb(i) = exp(-r*T)*mean(max(S-K,0));
    VpMCb(i) = exp(-r*T)*mean(max(K-S,0));

end
VcMC2B = mean(VcMCb);
VpMC2B = mean(VpMCb);

cputime_MC2B = toc;
fprintf('%20s%14.10f%14.10f%14.10f\n','Monte Carlo 2B',VcMC2B,VpMC2B,cputime_MC2B)


end


