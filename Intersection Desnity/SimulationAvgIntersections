clc; 
clear all;
%close all;

R = 50;
nB = 10;
R_0_vec = 20:5:25;


iterations = 1e4;
r_vec = R*(rand(iterations,nB));
theta_vec = 2*pi*rand(iterations,nB);

cc = [];

for i=1:length(R_0_vec)
    R_0 = R_0_vec(i);
    l_limit = R_0;
    out = [];
    for iter = 1:iterations
        ind = find(r_vec(iter,:)<=R_0);
        x_1 = -l_limit;
        y_1 = -cot(theta_vec(iter,ind))*x_1 + r_vec(iter,ind).*csc(theta_vec(iter,ind));
        x_2 = l_limit;
        y_2 = -cot(theta_vec(iter,ind))*x_2 + r_vec(iter,ind).*csc(theta_vec(iter,ind));
        x2 = [x_1 x_2];
        y2 = [y_1 y_2];
        y2 = reshape(y2,[length(ind),2]);
        intx = [];
        inty = [];
        for ii=1:size(y2,1)
            for jj=ii+1:size(y2,1)
                [qq, ll] = intersections(x2,y2(ii,:),x2,y2(jj,:));
                intx = [intx qq];
                inty = [inty ll];
            end
        end
        int = [intx; inty];
        distance = cell2mat(arrayfun(@(i) pdist([0,0;int(:,i)']),1:length(intx),UniformOutput=false));
        if sum(distance<=R_0)~=0
            out = [out sum(distance<=R_0)];
        else
            out = [out 0];
        end
    end
    cc = [cc mean(out)];
end
% mean(out)
