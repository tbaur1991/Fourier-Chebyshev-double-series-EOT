% prepare 3D plot
figure(1)
tiledlayout(1,1)
nexttile
ax = gca;
set(gcf,'Color','w')
% plot measurements
p1 = plot3(meas(1,:),meas(2,:),meas(3,:),'or','MarkerFaceColor','r');
axis equal; hold on
% plot reference
corner_pts = [-a_ref -b_ref 0; a_ref -b_ref 0; a_ref -b_ref h_ref; -a_ref -b_ref h_ref;...
              -a_ref b_ref 0; a_ref b_ref 0; a_ref b_ref h_ref; -a_ref b_ref h_ref]';
corner_pts = R_ref*corner_pts + pos_ref;
corner_pts = corner_pts';
p2 = plot3([corner_pts(1,1) corner_pts(2,1)],[corner_pts(1,2) corner_pts(2,2)],[corner_pts(1,3) corner_pts(2,3)],'k','LineWidth',2);
hold on
plot3([corner_pts(2,1) corner_pts(3,1)],[corner_pts(2,2) corner_pts(3,2)],[corner_pts(2,3) corner_pts(3,3)],'k','LineWidth',2)
plot3([corner_pts(3,1) corner_pts(4,1)],[corner_pts(3,2) corner_pts(4,2)],[corner_pts(3,3) corner_pts(4,3)],'k','LineWidth',2)
plot3([corner_pts(4,1) corner_pts(1,1)],[corner_pts(4,2) corner_pts(1,2)],[corner_pts(4,3) corner_pts(1,3)],'k','LineWidth',2)
plot3([corner_pts(2,1) corner_pts(6,1)],[corner_pts(2,2) corner_pts(6,2)],[corner_pts(2,3) corner_pts(6,3)],'k','LineWidth',2)
plot3([corner_pts(6,1) corner_pts(5,1)],[corner_pts(6,2) corner_pts(5,2)],[corner_pts(6,3) corner_pts(5,3)],'k','LineWidth',2)
plot3([corner_pts(5,1) corner_pts(1,1)],[corner_pts(5,2) corner_pts(1,2)],[corner_pts(5,3) corner_pts(1,3)],'k','LineWidth',2)
plot3([corner_pts(3,1) corner_pts(7,1)],[corner_pts(3,2) corner_pts(7,2)],[corner_pts(3,3) corner_pts(7,3)],'k','LineWidth',2)
plot3([corner_pts(7,1) corner_pts(8,1)],[corner_pts(7,2) corner_pts(8,2)],[corner_pts(7,3) corner_pts(8,3)],'k','LineWidth',2)
plot3([corner_pts(8,1) corner_pts(4,1)],[corner_pts(8,2) corner_pts(4,2)],[corner_pts(8,3) corner_pts(4,3)],'k','LineWidth',2)
plot3([corner_pts(8,1) corner_pts(5,1)],[corner_pts(8,2) corner_pts(5,2)],[corner_pts(8,3) corner_pts(5,3)],'k','LineWidth',2)
plot3([corner_pts(7,1) corner_pts(6,1)],[corner_pts(7,2) corner_pts(6,2)],[corner_pts(7,3) corner_pts(6,3)],'k','LineWidth',2)

% plot estimated cylindric shape
syms theta u
if strcmp(filter,'GAM')
    pos = X(1:3,k,j);
    or = X(4,k,j);
    h = c1(X(5,k,j),0,'lower');
    shape = X(6:end,k,j);
elseif strcmp(filter,'ERHM')
    pos = [X(1:2,k,j);X_line(1,k,j)];
    or = X(3,k,j);
    h = c1(X_line(2,k,j),0,'lower');
    shape = X(4:end,k,j);
end
r = fourier_chebychev_series(shape,theta,u,nu,nth);
x = r*cos(theta);
y = r*sin(theta);
z = u*h/2;
p3 = fsurf(x,y,z,[0 2*pi -1 1],'FaceAlpha',0.6,'FaceColor','#4DBEEE');
% translate surf plot
t = hgtransform('Parent',ax);
set(p3,'Parent',t);
m = makehgtform('translate',pos);
set(t,'Matrix',m)
% rotate surf plot
t = hgtransform('Parent',t);
set(p3,'Parent',t);
Rz = makehgtform('zrotate',or);
set(t,'Matrix',Rz);

% legend
leg = legend([p1,p2,p3],{'Measurement','Reference','Estimate'});
leg.Interpreter = "latex";
leg.Orientation = "horizontal";
leg.NumColumns = 3;
leg.Box = "off";
leg.Layout.Tile = "south";
leg.FontSize = 12;
% draw
drawnow
hold off