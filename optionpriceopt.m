%{
File to run a numerical estimation of optimal parameters for a given vector
of option prices
%}
%{
clear all
close all
clc
%}
%{
Inputs:
FVR_temp: vector of parameters in the following order:
r
q
lambda
eta
u
u0
rho
mu
%}

function S_MSE=optionpriceopt(FVR_temp, S_Struct)
%% find c hat, the option price according to the COS-FFT
N=1000; %precision of calculations

%set of parameters from the paper
r=0;
q=0;
lambda=1.5768;
eta = 0.5751;
u = 0.0398;
u0 = 0.0175;
rho = -0.5711;
mu=0; %the mean of the distribution

% set of option characteristics
S0=100; % price of the underlying asset
vectorS0 = 1:10:301; % a range for S0
K=100; % strike price
vectorK=1:10:301; % a range for K
call=0; % the input to test a call option
put=1; % the input if we want a put option
t0=0; %initial date
T=1; %maturity

%% c_hat for different parameters

FVr_temp(1)=r;
FVr_temp(2)=q;
FVr_temp(3)=lambda;
FVr_temp(4)=eta;
FVr_temp(5)=u;
FVr_temp(6)=u0;
FVr_temp(7)=rho;
FVr_temp(8)=mu;

S_MSE.FVr_oa(1) = c_hat(N,vectorS0,K,r,lambda,eta,u,u0,rho,T,t0,mu,put);
S_MSE.I_nc      = 0;%no constraints
S_MSE.FVr_ca    = 0;%no constraint array
S_MSE.I_no      = 1;%number of objectives (costs)

%{
1) calculate reference option prices - use COS-FFT
2) write/ include a function that calculates option prices for parameters
--> rewrite COS-FFT sucht that the input are theta and X_i
3) write objective function: (c-c^/hat)^2
4) implement differential evolution
%}
