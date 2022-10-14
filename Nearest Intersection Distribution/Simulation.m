clc; clear all;
%% Initilization
R = 50; % Radius of circle where lines are generated
nB = 10;    % Number of Lines

xt = 0; yt = 0;   % Coordinates of test point

L = 5000;   % Length of each line
iterations = 1e6;

t_vec = linspace(eps,200,400);

%% Nearest Intersection Distribution Calculation
ni_dist = [];
for nB = nB_vec
    prob_a = zeros(1,length(t_vec));
    prob = zeros(1,length(t_vec));
    for t_i = 1 : length(t_vec)
        t = t_vec(t_i);
        r_vec = R * (rand(iterations,nB));
        theta_vec = 2*pi *rand(iterations,nB);
        d_vec = (r_vec./cos(theta_vec)) + (rt./tan(pi/2 - theta_vec));
        for iter = 1:iterations
            distance = d_vec(iter,:);
            tt = distance >= xt-t & distance <=xt+t;
            if sum(tt==1)>=1
                prob(t_i) = prob(t_i)+ 1/iterations;
            end
        end
    end
    ni_dist = [ni_dist; prob];
end
