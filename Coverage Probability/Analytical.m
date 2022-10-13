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

theta_vec = linspace(0,2*pi,101);
d1_vec = eps:0.1:100;%sqrt((x_t+R)^2+R^2);
%% PGFL and CDF Calculation
for i = 1:length(x_t_vec)
    x_t = x_t_vec(i);
    r_vec = 0:1:R;
    PGFL_I_NI{i} = fnAvF(gamma,d1_vec,theta_vec,r_vec,x_t,lambda);

    cdf_nn = zeros(1,length(d1_vec));
    for j = 1:length(d1_vec)
        cdf_nn(j) = DistanceDist(d1_vec(j),x_t,R,lambda,r_vec,theta_vec);
    end
    cdf_all{i} = cdf_nn;
end
%% Coverage Probabailty Calculation
inner_term = exp(-gamma*N0./(K*P*(d1_vec.^(-alpha))));
for i = 1:length(x_t_vec)
    x_t = x_t_vec(i);
    pgfl = PGFL_I_NI{1,i};
    %% Area of Domain Band
    Area_D = arrayfun(@(x) Ad(x,x_t,R),d1_vec);
    Area_D(1) = eps;
    %% PDF of nearest neighbour
    pdf = 1-((cdf_all{1,i}./(2*pi*R)).^(nB));
    pdf_d1 = [eps, diff(pdf)./diff(d1_vec)];

    G_I_NI = ((pgfl(1,:) + pgfl(2,:))./(2*pi*R)).^(nB-1);
    G_I = ((pgfl(1,:)./(Area_D)));
    G_I(1)=1;

    cov(i) = trapz(d1_vec,pdf_d1.*inner_term.*G_I.*G_I_NI);
end