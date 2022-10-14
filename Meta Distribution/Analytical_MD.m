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
lambda_ap = 0.1;   % Intensity of PPP on lines
x_t_vec = 0:1:100;    % Test point
p_t_vec = 0:0.01:1; % Transmisson Probability

r_vec = 0:1:R;
theta_vec = linspace(0,2*pi,101);
d1_vec = eps:0.1:100;%sqrt((x_t+R)^2+R^2);

%% PGFL and CDF Calculation
PGFL_I_NI = [];
for i = 1:length(x_t_vec)
    x_t = x_t_vec(i);
    temp = [];
    for j = 1:length(p_t_vec)
        p_t = p_t_vec(j);
        %% Calculate M_{1} or M_{-1}
        pgfl = fnA_M1(gamma,d1_vec,theta_vec,r_vec,alpha,x_t,lambda_ap,p_t);    % M_{1}
%         pgfl = fnA_Delay(gamma,d1_vec,theta_vec,r_vec,alpha,x_t,lambda,p_t);    % M_{-1} 
        temp = [temp; pgfl];
    end
    PGFL_I_NI = [PGFL_I_NI; temp];

    cdf_nn = zeros(1,length(d1_vec));
    for j = 1:length(d1_vec)
        cdf_nn(j) = DistanceDist(d1_vec(j),x_t,lambda_ap,r_vec,theta_vec);
    end
    cdf_all{i} = cdf_nn;
end 

%% Meta Calcualtions
b=1; % For M_{1}, b=1 and M_{-1}, b=-1
noise_term = exp(-b*gamma*N0./(K*P*(d1_vec.^(-alpha))));
results = [];
for j=1:length(x_t_vec)
    x_t = x_t_vec(j);
    meta = zeros(1,length(p_t_vec));
    pgfl = PGFL_I_NI(2*(j-1)*length(p_t_vec)+1:2*j*length(p_t_vec),:);
    %% PDF of Nearest Neighbour
    pdf = 1-((cdf_all{1,j}./(2*pi*R)).^(nB));
    pdf_d1 = [eps, diff(pdf)./diff(d1_vec)];
    %% Area of Domain Band
    Area_D = arrayfun(@(x) Ad(x,x_t,R),d1_vec);
    Area_D(1) = eps;
    %% Meta Distribution
    for i = 1:length(p_t_vec)
        p_t = p_t_vec(i);
        pgfl_temp = pgfl((2*i)-1:2*i,:);
        pgfl_i_ni = ((pgfl_temp(1,:) + pgfl_temp(2,:))./(2*pi*R)).^(nB-1);
        pgfl_i = ((pgfl_temp(1,:)./(Area_D)));
        t4(1)=1;
        meta(i) = trapz(d1_vec,pgfl_i.*pgfl_i_ni.*pdf_d1.*noise_term);
    end
    results = [results; meta];
end

%% Succesfull Transmisson Density
i=1;    % choose i for x_t value i.e. i=1 means x_t = 0, i=51 means x_t=50
M_1 = results(i,:);
std = lambda_ap*M_1.*p_t_vec;

%% Mean Local Delay
% [~, ind] = arrayfun(@(i) min(results(i,1:end-1)./p_t_vec(1:end-1)),1:1:length(p_t_vec));
% D_p = p_t_vec(ind);
