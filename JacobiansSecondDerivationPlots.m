%%
clear
p = Geometry3D.SampleVisualDirections(100);

x = p(:,1);
y = p(:,2);
z = p(:,3);
zer = zeros(size(p(:,1)));

R(1,1,:) = x;
R(1,2,:) = -y;
R(1,3,:) = -z;
R(2,1,:) = y;
R(2,2,:) = 1 - ( y.^2 ) ./ (1 + x);
R(2,3,:) = -( y .* z ) ./ (1 + x) ;
R(3,1,:) = z ;
R(3,2,:) = -( y .* z ) ./ (1 + x) ;
R(3,3,:) = 1 -  ( z.^2 ) ./ (1 + x) ;

Jv(1,1,:) = x.^2 - 1;
Jv(1,2,:) = x.*y;
Jv(1,3,:) = x.*z;
Jv(2,1,:) = x.*y;
Jv(2,2,:) = y.^2 - 1;
Jv(2,3,:) = y.*z;
Jv(3,1,:) = x.*z;
Jv(3,2,:) = y.*z;
Jv(3,3,:) = z.^2 - 1;


Jw(1,1,:) = zer;
Jw(1,2,:) = z;
Jw(1,3,:) = -y;
Jw(2,1,:) = -z;
Jw(2,2,:) = zer;
Jw(2,3,:) = x;
Jw(3,1,:) = y ;
Jw(3,2,:) = -x;
Jw(3,3,:) = zer;


RJv(1,1,:) = zer;
RJv(1,2,:) = zer;
RJv(1,3,:) = zer;
RJv(2,1,:) = y;
RJv(2,2,:) = -1 +  ( y.^2 ) ./ (1 + x);
RJv(2,3,:) = ( y .* z ) ./ (1 + x) ;
RJv(3,1,:) = z ;
RJv(3,2,:) = ( y .* z ) ./ (1 + x) ;
RJv(3,3,:) = - 1 +  ( z.^2 ) ./ (1 + x) ;


RJw(1,1,:) = zer;
RJw(1,2,:) = zer;
RJw(1,3,:) = zer;
RJw(2,1,:) = -z;
RJw(2,2,:) = - ( y.*z ) ./ (1 + x);
RJw(2,3,:) = 1-(  z.^2 ) ./ (1 + x) ;
RJw(3,1,:) = y ;
RJw(3,2,:) = -1+( y .^2 ) ./ (1 + x) ;
RJw(3,3,:) = ( y.*z ) ./ (1 + x) ;

RRJv = pagemtimes(R,Jv);
RRJw = pagemtimes(R,Jw);


figure('color','w')
tiledlayout(1,3)
nexttile 
set(gca,'NextPlot','add')
plotSphere()
hl(1) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(R(1,1,:)), squeeze(R(2,1,:)), squeeze(R(3,1,:)), 'AutoScaleFactor', 0.5);
hl(2) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(R(1,2,:)), squeeze(R(2,2,:)), squeeze(R(3,2,:)), 'AutoScaleFactor', 0.5);
hl(3) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(R(1,3,:)), squeeze(R(2,3,:)), squeeze(R(3,3,:)), 'AutoScaleFactor', 0.5);
legend(hl,{'x','y','z'}, 'Location', 'northeast','box','off','fontsize',14);
nexttile
set(gca,'NextPlot','add')
plotSphere()
hl(1) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(Jv(1,1,:)), squeeze(Jv(1,2,:)), squeeze(Jv(1,3,:)), 'AutoScaleFactor', 0.5);
hl(2) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(Jv(2,1,:)), squeeze(Jv(2,2,:)), squeeze(Jv(2,3,:)), 'AutoScaleFactor', 0.5);
hl(3) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(Jv(3,1,:)), squeeze(Jv(3,2,:)), squeeze(Jv(3,3,:)), 'AutoScaleFactor', 0.5);
legend(hl,{'x','y','z'}, 'Location', 'northeast','box','off','fontsize',14);
nexttile
set(gca,'NextPlot','add')
plotSphere()
hl(1) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(Jw(1,1,:)), squeeze(Jw(1,2,:)), squeeze(Jw(1,3,:)), 'AutoScaleFactor', 0.5);
hl(2) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(Jw(2,1,:)), squeeze(Jw(2,2,:)), squeeze(Jw(2,3,:)), 'AutoScaleFactor', 0.5);
hl(3) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(Jw(3,1,:)), squeeze(Jw(3,2,:)), squeeze(Jw(3,3,:)), 'AutoScaleFactor', 0.5);

legend(hl,{'x','y','z'}, 'Location', 'northeast','box','off','fontsize',14);



figure('color','w')
tiledlayout(1,2)
nexttile
set(gca,'NextPlot','add')
plotSphere()

hl(1) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(RJv(2,1,:)), squeeze(RJv(2,2,:)), squeeze(RJv(2,3,:)), 'AutoScaleFactor', 0.5);
hl(2) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(RJv(3,1,:)), squeeze(RJv(3,2,:)), squeeze(RJv(3,3,:)), 'AutoScaleFactor', 0.5);
legend(hl,{'1','2'}, 'Location', 'northeast','box','off','fontsize',14);
nexttile
set(gca,'NextPlot','add')
plotSphere()

hl(1) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(RJw(2,1,:)), squeeze(RJw(2,2,:)), squeeze(RJw(2,3,:)), 'AutoScaleFactor', 0.5);
hl(2) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(RJw(3,1,:)), squeeze(RJw(3,2,:)), squeeze(RJw(3,3,:)), 'AutoScaleFactor', 0.5);

legend(hl,{'1','2'}, 'Location', 'northeast','box','off','fontsize',14);



%%
figure('color','w')
tiledlayout(2,3);
nexttile
set(gca,'NextPlot','add')
axis equal;
set(gca,'xlim',[-1 1], 'ylim',[-1 1])
quiver(y.*acosd(x)/90,z.*acosd(x)/90, squeeze(RJv(2,1,:)), squeeze(RJv(3,1,:)),'linewidth',2)
xlabel( 'u_1','fontsize',14);
ylabel( 'u_2','fontsize',14);
set(gca,'xtick',[],'ytick',[])

nexttile
set(gca,'NextPlot','add')
axis equal;
set(gca,'xlim',[-1 1], 'ylim',[-1 1])
quiver(y.*acosd(x)/90,z.*acosd(x)/90, squeeze(RJv(2,2,:)), squeeze(RJv(3,2,:)),'linewidth',2)
xlabel( 'u_1','fontsize',14);
ylabel( 'u_2','fontsize',14);
set(gca,'xtick',[],'ytick',[])

nexttile
set(gca,'NextPlot','add')
axis equal;
set(gca,'xlim',[-1 1], 'ylim',[-1 1])
quiver(y.*acosd(x)/90,z.*acosd(x)/90, squeeze(RJv(2,3,:)), squeeze(RJv(3,3,:)),'linewidth',2)
xlabel( 'u_1','fontsize',14);
ylabel( 'u_2','fontsize',14);
set(gca,'xtick',[],'ytick',[])


nexttile
set(gca,'NextPlot','add')
axis equal;
set(gca,'xlim',[-1 1], 'ylim',[-1 1])
quiver(y.*acosd(x)/90,z.*acosd(x)/90, squeeze(RJw(2,1,:)), squeeze(RJw(3,1,:)),'linewidth',2)
xlabel( 'u_1','fontsize',14);
ylabel( 'u_2','fontsize',14);
set(gca,'xtick',[],'ytick',[])

nexttile
set(gca,'NextPlot','add')
axis equal;
set(gca,'xlim',[-1 1], 'ylim',[-1 1])
quiver(y.*acosd(x)/90,z.*acosd(x)/90, squeeze(RJw(2,2,:)), squeeze(RJw(3,2,:)),'linewidth',2)
xlabel( 'u_1','fontsize',14);
ylabel( 'u_2','fontsize',14);
set(gca,'xtick',[],'ytick',[])

nexttile
set(gca,'NextPlot','add')
axis equal;
set(gca,'xlim',[-1 1], 'ylim',[-1 1])
quiver(y.*acosd(x)/90,z.*acosd(x)/90, squeeze(RJw(2,3,:)), squeeze(RJw(3,3,:)),'linewidth',2)
xlabel( 'u_1','fontsize',14);
ylabel( 'u_2','fontsize',14);
set(gca,'xtick',[],'ytick',[])



% hl = [];
% % hl(1) = quiver3(x,y,z, dazdx, dazdy, dazdz,'linewidth',2,'DisplayName', '\partial xyz / \partial u');
% % hl(2) = quiver3(x,y,z, deldx, deldy, deldz,'linewidth',2,'DisplayName', '\partial xyz / \partial v');
% % quiver3(x,y,z, x, y, z,'linewidth',1)
% 
% 
% 
% % quiver3(zer, zer, zer,p(:,1),p(:,2),p(:,3));
% hold
% %hl(1) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(J(1,1,:)), squeeze(J(1,2,:)), squeeze(J(1,3,:)), 'AutoScaleFactor', 0.5);
% %hl(2) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(J(2,1,:)), squeeze(J(2,2,:)), squeeze(J(2,3,:)), 'AutoScaleFactor', 0.5);
% %hl(3) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(J(3,1,:)), squeeze(J(3,2,:)), squeeze(J(3,3,:)), 'AutoScaleFactor', 0.5);
% 
% 
% %hl(4) = quiver3(p(:,1),p(:,2),p(:,3), -dazdx, -dazdy, -dazdz, 'AutoScaleFactor', 0.5);
% %hl(5) = quiver3(p(:,1),p(:,2),p(:,3), -deldx, -deldy, -deldz, 'AutoScaleFactor', 0.5);
% 
% hl(1) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(JJ(1,1,:)), squeeze(JJ(1,2,:)), squeeze(JJ(1,3,:)), 'AutoScaleFactor', 0.5);
% hl(2) = quiver3(p(:,1),p(:,2),p(:,3), squeeze(JJ(2,1,:)), squeeze(JJ(2,2,:)), squeeze(JJ(2,3,:)), 'AutoScaleFactor', 0.5);
% 
% hl(3) = quiver3(p(:,1),p(:,2),p(:,3), -squeeze(J2(1,1,:)), -squeeze(J2(1,2,:)), -squeeze(J2(1,3,:)), 'AutoScaleFactor', 0.5);
% hl(4) = quiver3(p(:,1),p(:,2),p(:,3), -squeeze(J2(2,1,:)), -squeeze(J2(2,2,:)), -squeeze(J2(2,3,:)), 'AutoScaleFactor', 0.5);
% legend(hl,'Location', 'northeast','box','off','fontsize',14);


function plotSphere()
R = 1; % radius of the eye
step = 1;
range = 90;
[az, el] = meshgrid(-2*range:step:2*range,-range:step:range); % azimuths and elevations to include
[xs,ys,zs] = Geometry3D.FickToSphere(deg2rad(az),deg2rad(el));

% subplot(2,5,[1 2 6 7],'nextplot','add')
% axis
% draw sphere
mesh(xs,ys,zs,'FaceAlpha', 0.5,'facecolor',0.8*[1 1 1],'EdgeColor','none');
xlabel(gca,'x')
ylabel(gca,'y')
zlabel(gca,'z')
% line of sight
line([1 R*1.5 ],[0 0 ],[0 0 ],'color','r','linewidth',2)

view(125,15);
axis equal;
end