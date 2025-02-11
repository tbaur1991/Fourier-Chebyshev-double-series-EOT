% *************************************************************************
%    I   SSSS   DDDD
%    I   S      D   D 
%    I   SSSS   D   D    Institut fr Systemdynamik
%    I      S   D   D
%    I   SSSS   DDDD   
% *************************************************************************
% Author: Tim Baur
% University: University of Applied Sciences "HTWG Konstanz"
% E-Mail: tbaur@htwg-konstanz.de
% Creation date: 06.02.2025
% *************************************************************************
% Content: In this simulation environment, the non-parametric 
% Fourer-Chebyshev 3D EOT filters are implemented in a static scenario. A 
% Monte Carlo simulation can be started by setting 'mc = true' and 'nMC' as 
% the number of runs to be conducted. By setting 'artificial_noise' as true 
% or false, an artificial measurement noise can be estimated for the 
% cylindrical shape representation to process measurements from the 
% interior of the object. Plots can be visualized by setting 'do_plot' as 
% true. By setting 'inside_outside_measurements' as true, measurements are
% generated from the boundary and the interior.
%
% *************************************************************************

close all
clear
clc

% simulation environment
nSamples = 100;
nMeas = 50;
mc = true;
do_plot = true;
inside_outside_measurements = false;

% filter settings. filter either 'GAM' or 'ERHM'
artificial_noise = false;
filter = 'ERHM';
n_upd = 20;

% triangular distributions
pd = makedist("Triangular","a",0,"b",1);

% set reference parameters
pos_ref = [2;1;0];
or_ref = 0;
R_ref = [cos(or_ref) -sin(or_ref) 0; sin(or_ref) cos(or_ref) 0; 0 0 1];
a_ref = 3;
b_ref = 1.5;
h_ref = 3;

% set number of coefficients for FCDS shape expansion
nu = 2;
nth = 7;
nCoeff = 1 + nu + nth + nu*nth;

% System state and evaluation allocation for cylindric shape estimation
if mc
    nMC = 100;
else
    nMC = 1;
end
if strcmp(filter,'GAM')
    nx = 5 + nCoeff;
elseif strcmp(filter,'ERHM')
    nx = 3 + nCoeff;
    % line system state
    X_line = zeros(2,nSamples,nMC);
    P_line = zeros(2,2,nSamples,nMC);
    Q_line = 1e-4*eye(2);
    % sigma point parameters line
    nx_line = 2 + n_upd;
    alpha = 0.1; beta = 2; kappa = 0;
    lambda = alpha^2*(nx_line + kappa) - nx_line;
    % calculate weights for sigma points
    wm(1) = lambda/(nx_line + lambda);
    wc(1) = lambda/(nx_line + lambda) + (1 - alpha^2 + beta);
    wm(2:2*nx_line + 1) = 1/(2*(nx_line + lambda));
    wc(2:2*nx_line + 1) = 1/(2*(nx_line + lambda));
else
    error('Wrong filter setting!')
end
X = zeros(nx,nSamples,nMC);
P = zeros(nx,nx,nSamples,nMC);
times = zeros(nSamples,nMC);
ious = zeros(nSamples,nMC);

% measurement and process noise
sig_r = 0.1;
Q = 1e-4*eye(nx);

%% start simulation
for j = 1:nMC
    j
    % artificial measurement noise 
    meanInX = 0; meanInY = 0; varInX = 0; varInY = 0;
    tao = 200;
    for k = 1:nSamples
        % generate new measurement
        thetas = 2*pi*rand(1,nMeas);
        us = rand(1,nMeas);
        ss = ones(1,nMeas); 
        if inside_outside_measurements
            ss(us > 0.7) = random(pd,1,numel(ss(us > 0.7)));
        end
        rs = zeros(1,nMeas);
        for i = 1:nMeas
            rs(i) = rec_contour(a_ref,b_ref,thetas(i));
        end
        meas = pos_ref + R_ref*[ss.*rs.*cos(thetas); ss.*rs.*sin(thetas); us*h_ref] + sig_r*randn(3,nMeas);
    
        if k == 1
            % do position initialization
            pos = zeros(3,1);
            pos(1:2) = mean(meas(1:2,:),2);
            pos(3) = min(meas(3,:));
            % init orientation
            c = cov(meas(1,:),meas(2,:));
            [v,d] = eig(c); [~,ii] = max(diag(d));
            or = atan2(v(2,ii),v(1,ii));
            % init extent
            R_rot = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1];
            mm = R_rot'*(meas - pos);
            b = abs(max(mm(2,:)));
            h = c1(max(mm(3,:)) - min(mm(3,:)),0,'lower');
            % do state initialization
            if strcmp(filter,'GAM')
                X(:,k,j) = [pos(1:2);pos(3)+h/2;or;h;4*b;zeros(nCoeff-1,1)];
                P(:,:,k,j) = blkdiag(c,eye(3),2*eye(nCoeff));
            elseif strcmpi(filter,'ERHM')
                X(:,k,j) = [pos(1:2);or;4*b;zeros(nCoeff-1,1)];
                P(:,:,k,j) = blkdiag(c,1,2*eye(nCoeff));
                X_line(:,k,j) = [pos(3)+h/2;h];
                P_line(:,:,k,j) = eye(2);
            end
        else
            % do random walk prediction
            X(:,k,j) = X(:,k-1,j);
            P(:,:,k,j) = P(:,:,k-1,j) + Q;
            if strcmp(filter,'ERHM')
                X_line(:,k,j) = X_line(:,k-1,j);
                P_line(:,:,k,j) = P_line(:,:,k-1,j) + Q_line;
            end
        end
    
        % do measurement updates
        if strcmp(filter,'GAM')
            tic
            [X(:,k,j),P(:,:,k,j),meanInX,meanInY,varInX,varInY] = update_FCDS_GAM(X(:,k,j),P(:,:,k,j),...
               sig_r,meas,n_upd,nu,nth,artificial_noise,meanInX,meanInY,varInX,varInY,tao);
            times(k,j) = toc;
        elseif strcmp(filter,'ERHM')
            tic
            [X_line(:,k,j),P_line(:,:,k,j)] = update_line_ERHM(X_line(:,k,j),P_line(:,:,k,j),sig_r,meas,n_upd,lambda,wm,wc);
            [X(:,k,j),P(:,:,k,j),meanInX,meanInY,varInX,varInY] = update_FCDS_ERHM(X(:,k,j),P(:,:,k,j),X_line(:,k,j),...
               sig_r,meas,n_upd,nu,nth,artificial_noise,meanInX,meanInY,varInX,varInY,tao);
            times(k,j) = toc;
        end
    
        % plot reference, measurement and estimate
        if do_plot
            plot_results
        end

        % calculate intersection over union
        if strcmp(filter,'GAM')
            ious(k,j) = iouFCDS(X(:,k,j),pos_ref,or_ref,a_ref,b_ref,h_ref,nu,nth,20);
        elseif strcmp(filter,'ERHM')
            XX = [X(1:2,k,j);X_line(1,k,j);X(3,k,j);X_line(2,k,j);X(4:end,k,j)];
            ious(k,j) = iouFCDS(XX,pos_ref,or_ref,a_ref,b_ref,h_ref,nu,nth,20);
        end
    end
end

% plot final results
figure(2)
tiledlayout(3,1)
set(gcf,'Color','w')

% mean values
X_mean = mean(X,3);
iou_mean = mean(ious,2);

if strcmp(filter,'GAM')
    % constrained values
    X(5,:,:) = c1(X(5,:,:),0,'lower');
    % calculate rmses
    rmse_or = sqrt((X_mean(4,:) - or_ref).^2);
    rmse_h = sqrt((X_mean(5,:) - h_ref).^2);
elseif strcmp(filter,'ERHM')
    % constrained values
    X_line(2,:,:) = c1(X_line(2,:,:),0,'lower');
    % mean values
    X_line_mean = mean(X_line,3);
    % calculate rmses
    rmse_or = sqrt((X_mean(3,:) - or_ref).^2);
    rmse_h = sqrt((X_line_mean(2,:) - h_ref).^2);
end

%% plot position rmse
nexttile
plot(1:nSamples,iou_mean,'LineWidth',2)
ylabel('IoU','Interpreter','latex')
title('Intersection over Union')

% axes
set(gca,'TickLength',[0 0])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12)

%% plot orientation rmse
nexttile
plot(1:nSamples,rmse_or,'LineWidth',2)
ylabel('RMSE','Interpreter','latex')
title('Orientation RMSE')

% axes
set(gca,'TickLength',[0 0])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12)

%% plot height rmse
nexttile
plot(1:nSamples,rmse_h,'LineWidth',2)
xlabel('sample','Interpreter','latex')
ylabel('RMSE','Interpreter','latex')
title('Height RMSE')

% axes
set(gca,'TickLength',[0 0])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12)

% display mean calculation time
disp(['Mean calculation time FCDS: ',num2str(mean(times(:))/nMeas)])

% end of function code
% *************************************************************************
%
%
% *************************************************************************