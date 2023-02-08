clc; clear all;

R=50;
nB=10;

r_vec = 0:0.5:100;

ind1 = r_vec<=R;
r2 = r_vec(~ind1);
lld_1 = (nB*(nB-1)/(4*pi*R^2))*ones(1,sum(ind1==1));
lld_2 = (nB*(nB-1)./(4*pi^2*R^2*r2)).*(2*r2.*asin(R./r2) - (2*R./r2).*sqrt(r2.^2 - R^2));
lld = [lld_1 lld_2];
