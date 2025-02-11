function r = rec_contour(a,b,theta)
    if theta < 0
        theta = 2*pi + theta;
    end
    if theta > 2*pi
        theta = theta - 2*pi;
    end
    phi = atan(b/a);
    if theta > 2*pi-phi || theta >=0 && theta <= phi
        r = a/cos(theta);
    elseif theta > phi && theta <= pi-phi
        r = b/sin(theta);
    elseif theta > pi-phi && theta <= pi+phi
        r = a/cos(pi-theta);
    elseif theta > pi+phi && theta <= 2*pi-phi
        r = b/sin(theta-pi);
    end
end

