function H = jacobi_meas_FCDS_ERHM(theta,u,r,pos,or,shape,meas,nu,nth)
    % get chebychev function values
    Ts = chebyshevList(nu,u);
    % derivative to position
    R = [cos(or) -sin(or); sin(or) cos(or)];
    z_loc = R'*(meas - pos);
    dr_dtheta = deriv_dr_dtheta(shape,theta,nu,nth,Ts);
    dr_du = deriv_dr_du(shape,theta,u,nu,nth);
    dtheta_dzloc = [-z_loc(2)/(z_loc(1)^2+z_loc(2)^2) z_loc(1)/(z_loc(1)^2+z_loc(2)^2)];
    dtheta_dm = -dtheta_dzloc*R';
    du_dm = -[0 0]*R';
    dr_dm = [dr_dtheta dr_du]*[dtheta_dm; du_dm];
    dcos_dm = -sin(theta)*dtheta_dm;
    dsin_dm = cos(theta)*dtheta_dm;
    dhx_dpos = eye(2) + R*[dr_dm*cos(theta) + r*dcos_dm; dr_dm*sin(theta) + r*dsin_dm];
    % derivative to orientation
    dR_dor = [-sin(or) -cos(or); cos(or) -sin(or)];
    dzloc_dor = dR_dor'*(meas - pos);
    dtheta_dor = dtheta_dzloc*dzloc_dor;
    du_dor = [0 0]*dzloc_dor;
    dr_dor = [dr_dtheta dr_du]*[dtheta_dor; du_dor];
    dcos_dor = -sin(theta)*dtheta_dor;
    dsin_dor = cos(theta)*dtheta_dor;
    dcart_dor = [dr_dor*cos(theta) + r*dcos_dor; dr_dor*sin(theta) + r*dsin_dor];
    dhx_dor = dR_dor*[r*cos(theta); r*sin(theta)] + R*dcart_dor;
    % derivative to shape parameters
    dr_dshape = deriv_dr_dshape(theta,nth,nu,Ts);
    dhx_dshape = R*[dr_dshape*cos(theta); dr_dshape*sin(theta)];
    % build jacobi matrix
    H = [dhx_dpos dhx_dor dhx_dshape];
end

function dr_dtheta = deriv_dr_dtheta(shape,theta,nu,nth,Ts)
    % get coefficients
    a0m = shape(nu+2:nu+nth+1);
    anm = shape(nu+nth+2:end);

    % second sum with vertical plane of symmetry
    s2 = -0.5*(1:nth).*sin((1:nth)*theta)*a0m;
    
    % third sum with vertical plane of symmetry
    s3 = 0;
    for m = 1:nth
        am = anm((m-1)*nu+1:m*nu);
        s3 = s3 - m*Ts*am*sin(m*theta);
    end
	
    % calculate radius 
    dr_dtheta = s2 + s3;
end

function dr_du = deriv_dr_du(shape,theta,u,nu,nth)
    % get coefficients
    an0 = shape(2:nu+1);
    anm = shape(nu+nth+2:end);

    % get chebychev derivative function values
    dTs = chebyshevListDeriv(nu,u);

    % first sum
    s1 = 0.5*dTs*an0;
    
    % third sum with vertical plane of symmetry
    s3 = 0;
    for m = 1:nth
        am = anm((m-1)*nu+1:m*nu);
        s3 = s3 + dTs*am*cos(m*theta);
    end
	
    % calculate radius 
    dr_du = s1 + s3;
end

function dr_dshape = deriv_dr_dshape(theta,nth,nu,Ts)    
    % start derivative array
    dr_dshape = [1/4 0.5*Ts];

    % derivative array with vertical plane of symmetry
    dr_dshape = [dr_dshape 0.5*cos((1:nth)*theta)];

    % derivative array with vertical plane of symmetry
    dr_dshape2 = zeros(1,nth*nu);
    for m = 1:nth
        dr_dshape2((m-1)*nu+1:m*nu) = Ts*cos(m*theta);
    end
    dr_dshape = [dr_dshape dr_dshape2];
end