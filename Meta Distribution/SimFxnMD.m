function out = SimFxnMD(gamma_vec,iterations,nB,xt,R,L,P,K,N0,alpha,lambda,p_t,rel_th)
    out = [];
    meta = zeros(1,iterations);
    for gamma_idx = 1:length(gamma_vec)
        r_vec = R * rand(iterations,nB);
        theta_vec = 2*pi*rand(iterations,nB);
        line_d_vec = abs(xt * cos(theta_vec) - r_vec);
        N_points = poissrnd(lambda*L,iterations, nB);
        for iter = 1:iterations
            point_dist = arrayfun(@(point) L*rand(1,point) - (L/2), N_points(iter,:), 'UniformOutput', 0);
            distance = cell2mat(arrayfun(@(lines)  sqrt((line_d_vec(iter,lines))^2 ...
                + (cell2mat(point_dist(lines))).^2) , 1:nB, 'un', 0));
            [d1,idx] = min(distance);
            interference = distance;
            interference(idx) = [];
            interference = (interference.^(-alpha)).*exprnd(1,size(interference));
            noise_term_suc_prob = exp(-gamma_vec(gamma_idx)*N0./(K*P*d1^(-alpha)));
            meta(iter) = noise_term_suc_prob*prod(p_t./(1 + gamma_vec(gamma_idx)*interference./(d1^(-alpha))) + (1-p_t));
        end

        for u = 1:length(rel_th)
            md_vs_rth(u) = mean(meta > rel_th(u));
        end
        out = [out; md_vs_rth];
    end
end