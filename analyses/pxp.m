%1. Load negative log likelihoods from a given experiment
m=csvread('modelResults/Exp1diffevidence.csv'); %Experiment 1
%m=csvread('modelResults/Exp2diffevidence.csv'); %Experiment 2

%Compute the values from the negative loss
[alpha,exp_r,xp,protectedExceedenceProb,bor] = bms(-m);
%probability of exceedance
probofexp=protectedExceedenceProb;
%Bayesian omnibus risk
bayesianomni=bor;
%Save results
csvwrite('modelResults/Exp1PXP.csv', protectedExceedenceProb) %Experiment 1
%csvwrite('modelResults/Exp2PXP.csv', protectedExceedenceProb) %Experiment 2