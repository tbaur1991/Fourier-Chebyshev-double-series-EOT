function r = fourier_chebychev_series(shape,theta,u,nu,nth)
    % get coefficients
    a00 = shape(1);
    an0 = shape(2:nu+1);
    a0m = shape(nu+2:nu+nth+1);
    anm = shape(nu+nth+2:end);

    % get chebychev function values
    Ts = chebyshevList(nu,u);

    % first sum
    s1 = 0.5*Ts*an0;

    % second sum with vertical plane of symmetry
    s2 = 0.5*cos((1:nth)*theta)*a0m;
    
    % third sum with vertical plane of symmetry
    s3 = 0;
    for m = 1:nth
        am = anm((m-1)*nu+1:m*nu);
        s3 = s3 + Ts*am*cos(m*theta);
    end
	
    % calculate radius 
    r = a00/4 + s1 + s2 + s3;
end