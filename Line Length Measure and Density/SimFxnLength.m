function out = SimFxnLength(gamma_vec,iterations,nB,xt,R,t)
    for gamma_idx = 1:length(gamma_vec)

        r_vec = R * rand(iterations,nB);
        theta_vec = 2*pi*rand(iterations,nB);
        line_d_vec = abs(xt*cos(theta_vec) - r_vec);
        ind = arrayfun(@(i) (line_d_vec(i,:)<=t),1:iterations,'UniformOutput',false);
        temp = [];
        for iter = 1:iterations
            ind2 = find(ind{iter}==1);
            if ~isempty(ind2)
                r = line_d_vec(iter,ind2);
                distance = 2*sqrt(t^2 - r.^2);
                temp = [temp sum(distance)];
            else
                temp = [temp 0];
            end
        end
    end
    out = temp;
end
