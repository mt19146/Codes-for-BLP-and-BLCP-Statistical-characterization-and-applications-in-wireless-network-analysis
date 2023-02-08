clc; clear all;

R=50;
nB=10;

r_vec = 0:0.5:100;

ind1 = r_vec<=R;
r2 = r_vec(~ind1);

lld_1 = (nB/(2*R))*ones(1,sum(ind1==1));
lld_1 = (nB/(pi*R)).*asin(R./r2);
lld = [lld_1 lld_2];
