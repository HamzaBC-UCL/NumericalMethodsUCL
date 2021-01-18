%% Instructions / ReadMe

% You can run all and outputs are printed in the command window
% Or you can also step through each section
% Each section calls a function to produce the output
% Each function & section contains some comments explaining my thought process.

%% Parameters for each function
clear

% Market parameters
T = 1; % maturity
S0 = 1; % spot price
K = 1.1; % strike price
r = 0.05; % risk-free interest rate
q = 0.02; % dividend rate

% Model parameter
sigma = 0.4; % volatility

% Risk-neutral measure
muRN = r-q-0.5*sigma^2; % drift

% Monte Carlo parameters; npaths = nblocks*nsample
nblocks = 20000; % number of blocks
nsample = 10000; % number of samples per block

% Fourier parameters
xwidth = 6; % width of the support in real space
ngrid = 2^8; % number of grid points
alpha = -10; % damping factor for a call

%% Q1

%Calling the BSM function

AnalyticalBSM(S0,K,T,sigma,q,r);

% Calling the Fourier transform function
fourierBSM(ngrid,xwidth,alpha,muRN,sigma,T,S0,K,r);

% Calling the Monte Carlo function
mcBSM(nblocks,muRN,r,sigma,nsample,T,S0,K);

%% Q2
alpha2A = -1;
muj =  -0.1;
sigmaj = 0.15;
lamda = 0.5;
muRN2A = r-q-0.5*sigma^2 - lamda*(exp(muj + 0.5*(sigmaj)^2)-1)  ; % RN drift from SDE

% 2A)Fourier with jumps
[VcF2A,VpF2A,cputime_F2A] = fourier2A(ngrid,xwidth,alpha2A,muRN2A,sigma,T,S0,K,r,muj,sigmaj,lamda);

%% 2B) Monte Carlo with jumps

[VcMC2B,VpMC2B,cputime_MC2B] = mc2B(nblocks,muRN2A,r,sigma,nsample,T,S0,K,lamda,muj,sigmaj);

%% Q3

% 3A, complex conjugate integrand

[VcF3A,VpF3A,cputime_F3A] = fourier3A(ngrid,xwidth,alpha2A,muRN2A,sigma,T,S0,K,r,muj,sigmaj,lamda);

%% 3B, trapezodial integration FFT for jumps
% I wasnt sure whether to use the normal GBM or jump, it works for both either way

[VcF3B,VpF3B,cputime_F3B] = fourier_trapz(ngrid,xwidth,alpha2A,muRN2A,sigma,T,S0,K,r,muj,sigmaj,lamda);

%% 3C, sum FFT

[VcF3C,VpF3C,cputime_F3C] = fourier_sum(ngrid,xwidth,alpha2A,muRN2A,sigma,T,S0,K,r,muj,sigmaj,lamda);


%% 3D PDFs, note that the MC PDF changes on each run because of the rand

% 3D i - Monte Carlo PDF hist, using MC with jumps
[X] = mcHistPDF(nblocks,muRN2A,sigma,nsample,T,S0,lamda,muj,sigmaj);

% 3D ii - inverse fourier PDF, using fourier with jumps
[pdf,x] = fourierPDF(ngrid,xwidth,muRN2A,sigma,T,muj,sigmaj,lamda);

% sclaing for axis', for line chart version of MC PDF
x_axis = linspace(max(x),min(x),ngrid/2);
bins = hist(X,ngrid/2);
mc_pdf = bins/sum(bins);

hold on
histogram(X,ngrid/2,'normalization','probability','FaceColor','g');

% For the line chart MC PDF, makes the integration check easier below
%plot(x_axis,mc_pdf,'g'); 

plot(x,pdf,'b','LineWidth',2);

title(sprintf('Monte Carlo and Inverse FFT PDFs'))
legend('Monte Carlo','Inverse FFT')
ylabel('Prob. %')
xlabel('Value')

hold off

% Checking if area underneath PDFs sum to 1
% fprintf('%20s%14.10f%14.10f%14.10f\n','PDF Integration Test FFT',trapz(pdf))
% fprintf('%20s%14.10f%14.10f%14.10f\n','PDF Integration Test MC',trapz(mc_pdf))
%% 3E

% I used the new payoff function, which removes the discounting factor
% Intuitively, This causes the prices to increase by a lot
% with the same payoff as previous questions the prices are similar to the
% other questions.

fourier3E(ngrid,xwidth,alpha2A,muRN2A,sigma,T,S0,K,r,muj,sigmaj,lamda);
