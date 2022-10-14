function out = fnA_M1(gamma_vec,d1_vec,theta_vec,r_vec,alpha,x_t,lambda,p_t)
    gamma_all = [];
    for gamma = gamma_vec
        d1_all1 = [];
        d1_all2 = [];
        for d1 = d1_vec
            for i = 1:length(theta_vec)
                theta = theta_vec(i);
                valid_C = find(d1 >= abs(x_t.*cos(theta) - r_vec));
                Nvalid_C = find(d1 < abs(x_t.*cos(theta) - r_vec));
                first_I1 = zeros(1,length(r_vec));
                first_I2 = zeros(1,length(r_vec));
                r_vec1 = r_vec(valid_C);
                r_vec2 = r_vec(Nvalid_C);

                first_I1(valid_C) = exp(-2*p_t*lambda*gamma*((d1^alpha)./(sqrt(gamma*(d1^alpha) + ...
                            (x_t*cos(theta) - r_vec1).^2))).*((pi/2) - atan(sqrt(((d1^alpha) - ...
                            (x_t*cos(theta) - r_vec1).^2)./(gamma*(d1^alpha) + (x_t*cos(theta) - r_vec1).^2)))));

                first_I1(Nvalid_C) = 0;

                first_I2(Nvalid_C) = exp(-2*p_t*lambda*pi*gamma*(d1^alpha)./(2*sqrt(gamma*(d1^alpha) + (x_t*cos(theta) - r_vec2).^2)));
                first_I2(valid_C) = 0;

                second_I1(i) = trapz(r_vec,first_I1);
                second_I2(i) = trapz(r_vec,first_I2);
            end
            theta_avg1 = trapz(theta_vec,second_I1);
            theta_avg2 = trapz(theta_vec,second_I2);
            d1_all1 = [d1_all1 theta_avg1];
            d1_all2 = [d1_all2 theta_avg2];
        end
        gamma_all = [gamma_all; [d1_all1; d1_all2]]; 
    end
    out = gamma_all;
end
