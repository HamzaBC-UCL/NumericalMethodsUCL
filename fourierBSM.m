function [VcF,VpF,cputime_F] = fourierBSM(ngrid,xwidth,alpha,muRN,sigma,T,S0,K,r)

% Function that produces the fourier prices of a call and put option
% + CPU time

% Grids in real and Fourier space
tic
N = ngrid/2;
b = xwidth/2; % upper bound of the support in real space
dx = xwidth/ngrid;
x = dx*(-N:N-1);
dxi = 2*pi/xwidth; % Nyquist relation
xi = dxi*(-N:N-1);

% Characteristic function at time T
xia = xi+1i*alpha; % call
psi = 1i*muRN*xia-0.5*(sigma*xia).^2; % characteristic exponent
Psic = exp(psi*T); % characteristic function
xia = xi-1i*alpha; % put
psi = 1i*muRN*xia-0.5*(sigma*xia).^2; % characteristic exponent
Psip = exp(psi*T); % characteristic function


% Fourier transform of the payoff
U = S0*exp(b);
L = S0*exp(-b);
[~,gc,Gc] = payoff(x,xi,alpha,K,L,U,S0,1); % call
[S,gp,Gp] = payoff(x,xi,-alpha,K,L,U,S0,0); % put

% Discounted expected payoff computed with the Plancherel theorem
c = exp(-r*T).*real(fftshift(fft(ifftshift(Gc.*conj(Psic)))))/xwidth; % call
VcF = interp1(S,c,S0,'spline');

p = exp(-r*T).*real(fftshift(fft(ifftshift(Gp.*conj(Psip)))))/xwidth; % put
VpF = interp1(S,p,S0,'spline');

cputime_F = toc;

fprintf('%20s%14.10f%14.10f%14.10f\n','Fourier',VcF,VpF,cputime_F)

end




