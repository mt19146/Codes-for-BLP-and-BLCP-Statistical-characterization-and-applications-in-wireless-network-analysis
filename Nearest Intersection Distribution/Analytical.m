clc; clear all;
%% Input
R=50;
nB=5;
r_0 = 100;

t_vec = linspace(eps,100,200);
theta = linspace(0,pi,500);
alpha_vec = linspace(0,pi,500);

prob = [];
tt=[];

for i = 1:length(t_vec)
    t=t_vec(i);
    prob1=[];
    for j=1:length(alpha_vec)
        alpha=alpha_vec(j);
        if alpha<=(pi/2)
            r_L = [max(-R,min(R,r_0*cos(theta(theta<=alpha+pi/2)) - t*cos(theta(theta<=alpha+pi/2)-alpha))) ...
                max(-R,min(R,r_0*cos(theta(theta>alpha+pi/2)) + t*cos(theta(theta>alpha+pi/2)-alpha)))];
            r_U = [max(-R,min(R,r_0*cos(theta(theta<=alpha+pi/2)) + t*cos(theta(theta<=alpha+pi/2)-alpha))) ...
                max(-R,min(R,r_0*cos(theta(theta>alpha+pi/2)) - t*cos(theta(theta>alpha+pi/2)-alpha)))];
        else
            r_L = [min(R,max(-R,r_0*cos(theta(theta<=alpha-pi/2)) + t*cos(theta(theta<=alpha-pi/2)-alpha))) ...
                max(-R,min(R,r_0*cos(theta(theta>alpha-pi/2)) - t*cos(theta(theta>alpha-pi/2)-alpha)))];
            r_U = [min(R,max(-R,r_0*cos(theta(theta<=alpha-pi/2)) - t*cos(theta(theta<=alpha-pi/2)-alpha))) ...
                max(-R,min(R,r_0*cos(theta(theta>alpha-pi/2)) + t*cos(theta(theta>alpha-pi/2)-alpha)))];
        end
        prob1 = [prob1 (trapz(theta,r_U) - trapz(theta,r_L))/(2*pi*R)];
    end
    prob = [prob trapz(alpha_vec, prob1)/(pi)];%[prob prob1];%
end
prob_o = 1-(1-prob./(1)).^nB;
