%{
File to run a numerical estimation of optimal parameters for a given vector
of option prices
%}

clear all
close all
clc

%% Input
S=50:10:150; % a vector of option prices, can be varied
%probably also strike price, time to maturity, price of underlying etc.?

%% Optimization

%{
we need: parameters that need to be determined --> mu, sigma, lambda and all that stuff?
         find out whether we need more input, e.g. a density/
            characteristic function
         the optimization algorithm --> find out how differential evolution
            is employed
%}
