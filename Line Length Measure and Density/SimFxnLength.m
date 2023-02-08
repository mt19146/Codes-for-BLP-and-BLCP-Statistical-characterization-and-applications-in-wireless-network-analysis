function out = SimFxnLength(gamma_vec,iterations,nB,xt,R,t)
    for gamma_idx = 1:length(gamma_vec)

        r_vec = R * rand(iterations,nB);
        theta_vec = 2*pi*rand(iterations,nB);
        line_d_vec = abs(xt*cos(theta_vec) - r_vec);
        ind1 = arrayfun(@(i) (line_d_vec(i,:)<=t(1)),1:iterations,'UniformOutput',false);
        ind2 = arrayfun(@(i) (line_d_vec(i,:)<=t(2)),1:iterations,'UniformOutput',false);
        temp = [];
        for iter = 1:iterations
            ind11 = find(ind1{iter}==1);
            ind22 = find(ind2{iter}==1);
            if ~isempty(ind11)
                r = line_d_vec(iter,ind11);
                distance = 2*sqrt(t(1)^2 - r.^2);
                temp1 = sum(distance);
            else
                temp1 = 0;
            end
            if ~isempty(ind22)
                r = line_d_vec(iter,ind22);
                distance = 2*sqrt(t(2)^2 - r.^2);
                temp2 = sum(distance);
            else
                temp2 = 0;
            end
            temp = [temp temp2-temp1];
        end
    end
    out = temp;
end
