clc; clear all;

R = 50;
nB_vec = 5;%[5 10];
iterations = 1e4;
r_0 = 100;
l_limit = 300;


t_vec = linspace(eps,100,200);

out = [];
for nB = nB_vec
    prob_a = zeros(1,length(t_vec));
    prob = zeros(1,length(t_vec));
    for t_i = 1 : length(t_vec)
        t = t_vec(t_i);
        r_vec = R * (rand(iterations,nB));
        theta_vec = 2*pi *rand(iterations,nB);

        x_1 = -l_limit;
        y_1 = -cot(theta_vec)*x_1 + r_vec.*csc(theta_vec);
        
        x_2 = l_limit;
        y_2 = -cot(theta_vec)*x_2 + r_vec.*csc(theta_vec);
        
        x2 = [x_1 x_2];
        y2 = [y_1 y_2];
        for iter = 1:iterations
            alpha=pi*rand;
            x_1 = -l_limit;
            y_1 = tan(alpha)*(x_1-r_0);
            x_2 = l_limit;
            y_2 = tan(alpha)*(x_2-r_0);
            x1 = [x_1 x_2];
            y1 = [y_1 y_2];
            y2_s = y2(iter,:);
            [xv, ~] = arrayfun(@(i) intersections(x1,y1,x2,y2_s([i i+nB])),1:nB,UniformOutput=false);
            tf = cellfun('isempty',xv);
            xv(tf) = {inf};
            distance = cell2mat(xv);
            if alpha<=pi/2
                tt = distance >= r_0-t*cos(alpha) & distance <=r_0+t*cos(alpha);
            else
                tt = distance >= r_0+t*cos(alpha) & distance <=r_0-t*cos(alpha);
            end
            if sum(tt==1)>=1
                prob(t_i) = prob(t_i)+ 1/iterations;
            end
        end
    end
    out = [out; prob];
end
