function [VcF3A,VpF3A,cputime_F3A] = fourier3A(ngrid,xwidth,alpha2A,muRN2A,sigma,T,S0,K,r,muj,sigmaj,lamda)

% Function that produces the fourier prices of a call and put option
% + CPU time, taking the complex conjugate

% Grids in real and Fourier space
tic
N = ngrid/2;
b = xwidth/2; % upper bound of the support in real space
dx = xwidth/ngrid;
x = dx*(-N:N-1);
dxi = 2*pi/xwidth; % Nyquist relation
xi = dxi*(-N:N-1);

% Characteristic function at time T
xia = xi+1i*alpha2A; % call
psi = 1i*muRN2A*xia-0.5*(sigma*xia).^2 + (lamda*(exp(1i*muj*xia-0.5*(sigmaj*xia).^2)-1)); % characteristic exponent
Psic = exp(psi*T); % characteristic function

xia = xi-1i*alpha2A; % put
psi = 1i*muRN2A*xia-0.5*(sigma*xia).^2 + (lamda*(exp(1i*muj*xia-0.5*(sigmaj*xia).^2)-1)); % characteristic exponent
Psip = exp(psi*T); % characteristic function


% Fourier transform of the payoff
U = S0*exp(b);
L = S0*exp(-b);
[~,gc,Gc] = payoff(x,xi,alpha2A,K,L,U,S0,1); % call
[S,gp,Gp] = payoff(x,xi,-alpha2A,K,L,U,S0,0); % put

% Discounted expected payoff computed with the Plancherel theorem
c = exp(-r*T).*real(fftshift(fft(ifftshift(conj(Gc).*(Psic)))))/xwidth; % call
VcF3A = (interp1(S,c,S0,'spline'));


p = exp(-r*T).*real(fftshift(fft(ifftshift(conj(Gp).*(Psip)))))/xwidth; % put
VpF3A = (interp1(S,p,S0,'spline'));


cputime_F3A = toc;


fprintf('%20s%14.10f%14.10f%14.10f\n','Fourier 3A',VcF3A,VpF3A,cputime_F3A)

end




