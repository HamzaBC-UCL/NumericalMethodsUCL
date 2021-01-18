function [VcF3E,VpF3E,cputime_F3E] = fourier3E(ngrid,xwidth,alpha2A,muRN2A,sigma,T,S0,K,r,muj,sigmaj,lamda)

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
dt = T/(N*2); % scaling factor for final call and put values

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

% Call and Put Prices.
% With the new payoff function, we arent discoutning / dampening, therefore the
% put and option prices will be higher.

% Call
c = fftshift(fft(ifftshift((Psic)))); % get inverse FFT PDF for call
c = real(c)*dt; % scale by dt 

% discounting + integrating over the inversed PDF from FFT using trapz
VcF3E = trapz(max((S-K),0).*(c))*exp(-r*T); % use payoff from question

% Put
p = fftshift(fft(ifftshift((Psip)))); % get inverse FFT PDF for put
p = real(p)*dt; % scale by dt

% discounting + integrating over the inversed PDF from FFT usinf trapz
VpF3E = trapz(max((K-S),0).*(p))*exp(-r*T); % use payoff from question
cputime_F3E = toc;


fprintf('%20s%14.10f%14.10f%14.10f\n','Fourier 3E',VcF3E,VpF3E,cputime_F3E)
end




