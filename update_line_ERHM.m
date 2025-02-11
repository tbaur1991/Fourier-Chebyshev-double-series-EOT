function [X,P] = update_line_ERHM(X,P,sig_r,Y,n_upd,lambda,wm,wc)
    % get state dimension
    nx = length(X); nMeas = size(Y,2);
    nxu = nx + n_upd;

    % check if number of measurements smaller than n_upd
    alpha = 0.1; beta = 2; kappa = 0;
    if nMeas < n_upd
        n_upd = nMeas;
        nxu = nx + n_upd;
        lambda = alpha^2*(nxu + kappa) - nxu;
        wm = zeros(1,2*nxu+1); wc = wm;
        wm(1) = lambda/(nxu + lambda);
        wc(1) = lambda/(nxu + lambda) + (1 - alpha^2 + beta);
        wm(2:2*nxu + 1) = 1/(2*(nxu + lambda));
        wc(2:2*nxu + 1) = 1/(2*(nxu + lambda));
    end

    % iterate n times over n_upd measurements
    for j = 1:ceil(nMeas/n_upd)
        % get new measurement set
        if j <= floor(nMeas/n_upd)
            meas = Y(:,n_upd*j-n_upd+1:n_upd*j);
        else
            meas = Y(:,n_upd*(j-1)+1:end);
            % change number of measurements to be updated and UKF weights
            n_upd = mod(nMeas,n_upd);
            nxu = nx + n_upd;
            lambda = alpha^2*(nxu + kappa) - nxu;
            wm = zeros(1,2*nxu+1); wc = wm;
            wm(1) = lambda/(nxu + lambda);
            wc(1) = lambda/(nxu + lambda) + (1 - alpha^2 + beta);
            wm(2:2*nxu + 1) = 1/(2*(nxu + lambda));
            wc(2:2*nxu + 1) = 1/(2*(nxu + lambda));
        end
    
        % calculate sigma points
        Xu = [X; zeros(n_upd,1)];
        L = blkdiag(chol(P)',eye(n_upd));
        A = sqrt(nxu + lambda) * L;
        sig_points = [zeros(size(Xu)) -A A];
        sig_points = sig_points + repmat(Xu, 1, size(sig_points, 2));
        sig_points(3:end,:) = 0.5*(1 + erf(sig_points(3:end,:)/sqrt(2)));
        nSamples = size(sig_points,2);
    
        z_predict = zeros(n_upd,nSamples);
        for i = 1:nSamples
            % get sample quantities
            z = sig_points(1,i);
            h = c1(sig_points(2,i),0,'lower');
            us = sig_points(3:end,i);
            z_predict(:,i) = abs(meas(3,:)' - z) - us*h/2;
        end
    
        % get predicted measurements
        z_pred = sum(wm.*z_predict,2);
    
        % get update matrices
        S = zeros(n_upd);
        Psi = zeros(2,n_upd);
        for i = 1:size(sig_points,2)
            S = S + wc(i)*(z_predict(:,i) - z_pred)*(z_predict(:,i) - z_pred)';
            Psi = Psi + wc(i)*(sig_points(1:2,i) - X)*(z_predict(:,i) - z_pred)';
        end
        S = S + sig_r^2*eye(n_upd);
    
        % do UKF measurement update
        K = Psi/S;
        X = X + K*(zeros(n_upd,1) - z_pred);
        P = P - K*Psi';
        P = 0.5*(P+P');
    end
end