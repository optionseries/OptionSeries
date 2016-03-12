load('Data.mat')  % Load Data
%%
S=1;         % Stock price
K=1;         % Strike
t=252/252;         % Time in years
r=0.04;      % Interest rate
theta=0.04;  % Long-term variance
v=0.05;      % Initial variance
rho=-0.8;    % Correlation between two Brownian motions
kappa=6;     % Mean-reversion rate
eta=0.2;     % Volatility of volatility
M=0.001;     % Scale factor, default value 0.001

[vsum, vmat]=HestonSeries(u,S,K,t,r,theta,v,rho,kappa,eta,M);
vmat         % Matrix form of u_{ij} (eta/(1+eta))^i ((v-theta)/(1+v-theta))^j, for convergence confirmation
vsum         % Sum of vmat