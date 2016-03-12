load('Data.mat')   % Load data
%%
S=100;           % Stock price
K=100;           % Strike
t=1;             % Time in years
r=0.1;           % Interest rate
v=0.09;          % Variance
N=11;            % Order of approximation, maximum 11

AmericanSeries(d,u,S,K,t,r,v,N)