clc; clear all;

R=50;
nB=10;
R_0_vec = 1:1:75;

ind1 = R_0_vec<=R;

p=R_0_vec(ind1)/R;
q=1-p;
qq=1./(1-p);

id_1 = 0.25*nB*(p.^2)*(nB-1);

r=R_0_vec(~ind1);
id_1 = 0.5*2*(r.^2.*asin(R./r) + R*(2*R*acos(R./r) - sqrt(r.^2 - R^2)))/(2*pi*R^2)*nB*(nB-1);

id = [id_1 id_2];
