clc; clear all;
%% Initilization
R = 50; % Radius of circle where lines are generated
nB = 10;    % Number of Lines

xt = 0; yt = 0;   % Coordinates of test point

t_vec = linspace(eps,200,400);
theta = linspace(eps,pi,500);

%%  Area of Domain bands corresponding to Nearest Intersection
ni_dist = [];

for i=1:length(t_vec)
    t=t_vec(i);

    ff1 = @(xt,t,yt,theta) sqrt(yt^2 + (xt-t)^2)*cos(atan(yt/(xt-t)) + theta);
    ff2 = @(xt,t,yt,theta) sqrt(yt^2 + (xt+t)^2)*cos(atan(yt/(xt+t)) + theta);
    
    r_L = ff1(xt,t,yt,theta);
    r_U = ff2(xt,t,yt,theta);
    
    
    r_L(r_L>R)=R;
    r_U(r_U>R)=R;
    r_L(r_L<-R)=-R;
    r_U(r_U<-R)=-R;
    
    
    if xt-t<=0
        r_L=-r_L;
    end
    if xt+t<=0
        r_U=-r_U;
    end
    
    temp = r_L;
    r_L(r_L>r_U) = r_U(r_L>r_U);
    r_U(temp>r_U) = temp(temp>r_U);
    ni_dist = [ni_dist (trapz(theta,r_U) - trapz(theta,r_L))];
end

%% Nearest Intersection Distribution Calculation
ni_dist = 1-(1-ni_dist./(2*R*pi)).^nB;
