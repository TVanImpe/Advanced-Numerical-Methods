%{
File to run a numerical estimation of optimal parameters for a given vector
of option prices
%}

clear all
close all
clc

%% Input
S=50:10:150; % a vector of option prices, can be varied

%% Optimization

%{
we need: parameters that need to be determined
         find out whether we need more input, e.g. a density/
            characteristic function
         the optimization algorithm --> find out how differential evolution
            is employed
%}
