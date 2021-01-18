function [Vca , Vpa,cputime_a] = AnalyticalBSM(S0,K,T,sigma,q,r)

% Function that produces the analytical prices of a BSM call and put option
% + CPU time

tic
d1 = (log(S0/K)+(r-q+0.5*sigma^2)*T)/(sigma*sqrt(T));
d2 = (log(S0/K)+(r-q-0.5*sigma^2)*T)/(sigma*sqrt(T));

Vca = S0*exp(-q*T)*cdf('Normal',d1,0,1) - K*exp(-r*T)*cdf('Normal',d2,0,1);
Vpa = K*exp(-r*T)*cdf('Normal',-d2,0,1) - S0*exp(-q*T)*cdf('Normal',-d1,0,1);

cputime_a = toc;


% Print the results
fprintf('%20s%14s%14s%14s\n','','call','put','CPU_time/s')
fprintf('%20s%14.10f%14.10f%14.10f\n','BS Analytical',Vca,Vpa,cputime_a)


end

