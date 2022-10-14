function P_cov = SimFxn(gamma_vec,iterations,nB,xt,R,L,P,K,N0,alpha,lambda)
    P_cov = zeros(1,length(gamma_vec));
    SINR = zeros(1,iterations);
    for gamma_idx = 1:length(gamma_vec)
        r_vec = R * rand(iterations,nB);
        theta_vec = 2*pi*rand(iterations,nB);
        line_d_vec = abs(xt*cos(theta_vec) - r_vec);
        N_points = poissrnd(lambda*L,iterations, nB);
        for iter = 1:iterations
            point_dist = arrayfun(@(point) L*rand(1,point) - (L/2), N_points(iter,:), 'UniformOutput', 0);
            distance = cell2mat(arrayfun(@(lines)  sqrt((line_d_vec(iter,lines))^2 ...
                + (cell2mat(point_dist(lines))).^2) , 1:nB, 'un', 0));
            [d1,idx] = min(distance);
            interference = distance;
            interference(idx) = [];
            interference2 = (interference.^(-alpha)).*exprnd(1,size(interference));
            x_d = (d1^-alpha)*exprnd(1);
            SINR(iter) = (P*K*x_d)/(N0 + (P*K*sum(interference2)));
        end
        P_cov(gamma_idx) = mean(SINR > gamma_vec(gamma_idx));
    end
end
