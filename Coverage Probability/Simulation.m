clc; clear all;
%% Initilization
P = 1;  % Power of Tx APs
alpha = 2;  % Fading Coeff.
B = 75e6;   % Bandwidth
N0 = B*10.^((-20.4));   % Noise
K = 1e-5;   % Path loss Constant
gamma = db2pow(-10); % Gamma threshold

%% Parameters
nB = 10;    % Number of Lines
R = 50; % Radius of circle where lines are generated
lambda = 0.1;   % Intensity of PPP on lines
x_t_vec = 0:1:2;    % Test point

L = 5000;   % Length of each line
iterations = 1e6;

%% Coverage Probabailty Calculation
for i=1:length(x_t_vec)
    xt = x_t_vec(i);
    P_cov(i) = SimFxn(gamma,iterations,nB,xt,R,L,P,K,N0,alpha,lambda);
end
