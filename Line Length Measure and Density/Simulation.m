clc; clear all;
L = 1000;

nB = 5;
R = 50;  
r0 = 0;

lambda = 0.1;

iterations = 1e5;

P = 1;
alpha = 2;

B = 75e6;
N0 = B*10.^((-20.4));
K = 1e-5;
gamma_vec = 0.1;

t_vec=1:1:50;
results = [];
for i=1:length(t_vec)
    t = t_vec(i);
    temp = SimFxnLength(gamma_vec,iterations,nB,r0,R,t);
    results = [results; temp];
end
