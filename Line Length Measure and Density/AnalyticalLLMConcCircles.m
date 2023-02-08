clc; clear all;
R=50;
xt=0;
nB=10;

t=1:1:75;
area_an = cell2mat(arrayfun(@(i) pi*(t(i+1))^2 - pi*(t(i))^2,1:length(t)-1,'UniformOutput',false));
L_bar = [pi*(t(t<=R))/2 sqrt(t(t>R).^2 - R^2) + asin(R./t(t>R)).*(t(t>R).^2)/R];

p = arrayfun(@(x) Ad(x,xt,R),t);
avg_length = diff(L_bar*nB.*(p/(2*pi*R)));

llm = avg_length./area_an;
