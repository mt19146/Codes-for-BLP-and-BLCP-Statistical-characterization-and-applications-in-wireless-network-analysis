clc; clear all;

R = 50;
R_0_vec = 1:1:75;

lambda = 0.1;
iterations = 1e5;

cc = [];

for i=1:length(R_0_vec)
    R_0 = R_0_vec(i);
    nB = poissrnd(lambda*2*pi*R_0,1,iterations);
    r_vec = arrayfun(@(i) R_0 * (rand(1,nB(i))),1:iterations,UniformOutput=false);
    theta_vec = arrayfun(@(i) 2*pi *rand(1,nB(i)),1:iterations,UniformOutput=false);
    l_limit = R_0+5;
    out = [];
    for iter = 1:iterations
        x_1 = -l_limit;
        y_1 = -cot(theta_vec{1,iter})*x_1 + r_vec{1,iter}.*csc(theta_vec{1,iter});
        x_2 = l_limit;
        y_2 = -cot(theta_vec{1,iter})*x_2 + r_vec{1,iter}.*csc(theta_vec{1,iter});
        x2 = [x_1 x_2];
        y2 = [y_1; y_2]';
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
