function [pdf,x] = fourierPDF(ngrid,xwidth,muRN2A,sigma,T,muj,sigmaj,lamda)

% Function that produces the fourier prices of a call and put option
% + CPU time + PDF

% Grids in real and Fourier space
N = ngrid/2;
dx = xwidth/ngrid;
x = dx*(-N:N-1);
dxi = 2*pi/xwidth; % Nyquist relation
xi = dxi*(-N:N-1);
dt = T/(N*2);

% Adjusted the char func such that we get the same PDF for either contract.
%
%%% Characteristic function at time T
psi = 1i*muRN2A*xi-0.5*(sigma*xi).^2 + (lamda*(exp(1i*muj*xi-0.5*(sigmaj*xi).^2)-1)); % characteristic exponent
Psic = exp(psi*T); % characteristic function

% PDF of inverse FFT characteristic function
pdf = real(fftshift(fft(ifftshift(Psic))))*dt;



end

