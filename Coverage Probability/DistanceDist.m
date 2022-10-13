function prob = DistanceDist(t,xt,lambda,r_vec,theta_vec)
    r_avg = zeros(size(r_vec));
    for r_i = 1:length(theta_vec)
        
        theta = theta_vec(r_i);
        C = zeros(size(r_vec));
        
        valid_chords = t >= abs(xt*cos(theta) - r_vec);
        C(valid_chords) = 2*sqrt(t^2 - (xt*cos(theta) - r_vec(valid_chords)).^2);
        void_prob = exp(-lambda*C);
        
        r_avg(r_i) = trapz(r_vec, void_prob);
    end
    prob = trapz(theta_vec,r_avg);
end


