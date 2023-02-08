function area = Ad(t,xt,R)
%Calculate Area of Domain Band

    if abs(xt) + t <= R % BLP = PLP
        area = 2*pi*t;
        
    elseif abs(xt) - t > R % Both Bands Clipped
        neg_term = (R-t)/xt;
        pos_term = (R+t)/xt;
        clipped_area = xt*(sqrt(1 - neg_term^2) - sqrt(1 - pos_term^2)) - (R-t)*acos(neg_term) + (R+t)*acos(pos_term);
        area = 2*pi*t - 2*clipped_area;
        
    else % One Band Clipped at a time
        neg_term = (R-t)/xt;
        clipped_area = xt*sqrt(1 - neg_term^2) - (R-t)*acos(neg_term);
        area = 2*pi*t - 2*clipped_area; 
    end
    if imag(area)~=0
        area = 2*pi*R;
    end
end

