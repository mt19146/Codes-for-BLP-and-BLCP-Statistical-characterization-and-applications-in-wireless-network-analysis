clc; clear all;
%% Initilization
P = 1;  % Power of Tx APs
alpha = 2;  % Fading Coeff.
B = 75e6;   % Bandwidth
N0 = B*10.^((-20.4));   % Noise
K = 1e-5;   % Path loss Constant
gamma_vec = db2pow([-20 -10 0]); % Gamma threshold

%% Parameters
nB = 10;    % Number of Lines
R = 50; % Radius of circle where lines are generated
lambda = 0.1;   % Intensity of PPP on lines
x_t_vec = 0;    % Test point
p_t= 0.99; % Transmisson Probability
rel_th = 0:0.01:1;  % Reliablity Threshold

L = 5000;   % Length of each line
iterations = 1e6;

%% Meta Distribution
results = [];
for i = 1:length(x_t_vec)
    xt = x_t_vec(i);
    results = [results; SimFxnMD(gamma_vec,iterations,nB,xt,R,L,P,K,N0,alpha,lambda,p_t,rel_th)];
end
