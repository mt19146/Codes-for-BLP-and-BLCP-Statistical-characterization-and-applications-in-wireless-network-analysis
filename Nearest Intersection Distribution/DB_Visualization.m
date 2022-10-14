clc; clear all;
%% Initilization
R = 50; % Radius of circle where lines are generated
nB = 10;    % Number of Lines

xt = 0; yt = 0;   % Coordinates of test point
t = 20;    %distance to nearest intersection

theta = linspace(eps,2*pi,500);

ff1 = @(xt,t,yt,theta) sqrt(yt^2 + (xt-t)^2)*cos(atan(yt/(xt-t)) + theta);
ff2 = @(xt,t,yt,theta) sqrt(yt^2 + (xt+t)^2)*cos(atan(yt/(xt+t)) + theta);

xx1 = ff1(xt,t,yt,theta);
xx2 = ff2(xt,t,yt,theta);

xx1(xx1>R)=R;
xx2(xx2>R)=R;
xx1(xx1<-R)=-R;
xx2(xx2<-R)=-R;   

if xt-t<=0
    xx1=-xx1;
end
if xt+t<=0
    xx2=-xx2;
end

figure;
plot(theta,xx1)
hold on
plot(theta,xx2)
hold off 
str = sprintf('R=%d, x_t=%d, t=%d, r_t=%d', R, xt, t, yt);
title(str)
xlabel("\theta [degrees]")
ylabel("r_L and r_U")
legend("r_L","r_U")
