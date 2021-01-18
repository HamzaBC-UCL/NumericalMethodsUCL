function [X] = mcHistPDF(nblocks,muRN2A,sigma,nsample,T,S0,lamda,muj,sigmaj)

% Function that produces merton jump monte carlo prices of a call and put option
% + CPU time

tic;

for i = 1:nblocks

    % Arithmetic Brownian motion 
    ADx = muRN2A*T + sigma*randn(1,nsample)*sqrt(T);
    
    % increments of NCPP
    dN = poissrnd(lamda*T,[1,nsample]);
    
    JDx = muj*dN + sigmaj*sqrt(dN).*randn(1,nsample); 
    
    % sum of increments NCPP & ABM
    X = JDx + ADx; % using this for the hist
    
    % not using the code below, kept it to play around with
    % Geometric Brownian motion at time T
    S = S0*exp(X);
    

end




end