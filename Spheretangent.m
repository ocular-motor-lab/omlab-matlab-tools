%%
% 
R = 1; % radius of the eye
step = 1;
range = 90;
[az, el] = meshgrid(-2*range:step:2*range,-range:step:range); % azimuths and elevations to include

[xs,ys,zs] = Geometry3D.FickToSphere(deg2rad(az),deg2rad(el));
a = ParticleSampleSphere( 'N' ,200);
x = a(a(:,1)>0,1);
y = a(a(:,1)>0,2);
z = a(a(:,1)>0,3);


dazdx = -y;
dazdy = 1 -  ( y.^2 ) ./ (1 + x);
dazdz = -( y .* z ) ./ (1 + x) ;

deldx = -z ;
deldy = -( y .* z ) ./ (1 + x) ;
deldz = 1 -  ( z.^2 ) ./ (1 + x) ;


rdazwx = 0*dazdx + z.*dazdy - y.*dazdz;
rdazwy = -z.*dazdx - 0.*dazdy + x.*dazdz;
rdazwz = y.*dazdx - x.*dazdy + 0*dazdz;

rdelwx = 0*deldx + z.*deldy - y.*deldz ;
rdelwy = -z.*deldx + 0*deldy + x.*deldz;
rdelwz = y.*deldx - x.*deldy + 0.*deldz ;



rdazwx = z;
rdazwy = ( y .* z ) ./ (1 + x);
rdazwz =  -1 +  ( z.^2 ) ./ (1 + x);

rdelwx = -y;
rdelwy = 1 -  ( y.^2 ) ./ (1 + x);
rdelwz = -( y .* z ) ./ (1 + x) ;



figure('color','w')
subplot(2,5,[1 2 6 7],'nextplot','add')
axis
% draw sphere
mesh(xs,ys,zs,'FaceAlpha', 0.5,'facecolor',0.8*[1 1 1],'EdgeColor','none');
xlabel(gca,'x')
ylabel(gca,'y')
zlabel(gca,'z')
% line of sight
line([1 R*1.5 ],[0 0 ],[0 0 ],'color','r','linewidth',2)

view(125,15);
axis equal;


hl = [];
hl(1) = quiver3(x,y,z, dazdx, dazdy, dazdz,'linewidth',2,'DisplayName', '\partial xyz / \partial u');
hl(2) = quiver3(x,y,z, deldx, deldy, deldz,'linewidth',2,'DisplayName', '\partial xyz / \partial v');
% quiver3(x,y,z, x, y, z,'linewidth',1)

legend(hl,'Location', 'northeast','box','off','fontsize',14);




% figure('color','w')
subplot(2,5,3,'nextplot','add')
axis equal;
set(gca,'xlim',[-1 1], 'ylim',[-1 1])
quiver(y.*acosd(x)/90,z.*acosd(x)/90, dazdx, deldx,'linewidth',2)
xlabel( '\partial u / \partial x','fontsize',14);
ylabel( '\partial v / \partial x','fontsize',14);
set(gca,'xtick',[],'ytick',[])

subplot(2,5,4,'nextplot','add')
axis equal;
set(gca,'xlim',[-1 1], 'ylim',[-1 1])
quiver(y.*acosd(x)/90,z.*acosd(x)/90, dazdy, deldy,'linewidth',2)
xlabel( '\partial u / \partial y','fontsize',14);
ylabel( '\partial v / \partial y','fontsize',14);
set(gca,'xtick',[],'ytick',[])

subplot(2,5,5,'nextplot','add')
axis equal;
set(gca,'xlim',[-1 1], 'ylim',[-1 1])
quiver(y.*acosd(x)/90,z.*acosd(x)/90, dazdz, deldz,'linewidth',2)
xlabel( '\partial u / \partial z','fontsize',14);
ylabel( '\partial v / \partial z','fontsize',14);
set(gca,'xtick',[],'ytick',[])


subplot(2,5,8,'nextplot','add')
axis equal;
set(gca,'xlim',[-1 1], 'ylim',[-1 1])
quiver(y.*acosd(x)/90,z.*acosd(x)/90, rdazwx, rdelwx,'linewidth',2)
xlabel( '\partial u / \omega_x dt','fontsize',14);
ylabel( '\partial v / \omega_x dt','fontsize',14);
set(gca,'xtick',[],'ytick',[])

subplot(2,5,9,'nextplot','add')
axis equal;
set(gca,'xlim',[-1 1], 'ylim',[-1 1])
quiver(y.*acosd(x)/90,z.*acosd(x)/90, rdazwy, rdelwy,'linewidth',2)
xlabel( '\partial u / \omega_y dt','fontsize',14);
ylabel( '\partial v / \omega_y dt','fontsize',14);
set(gca,'xtick',[],'ytick',[])

subplot(2,5,10,'nextplot','add')
axis equal;
set(gca,'xlim',[-1 1], 'ylim',[-1 1])
quiver(y.*acosd(x)/90,z.*acosd(x)/90, rdazwz, rdelwz,'linewidth',2)
xlabel( '\partial u / \omega _z dt','fontsize',14);
ylabel( '\partial v / \omega_z dt','fontsize',14);
set(gca,'xtick',[],'ytick',[])


%%


%%
% 
R = 1; % radius of the eye
step = 10;
range = 90;
[az, el] = meshgrid(-2*range:step:2*range,-range:step:range); % azimuths and elevations to include

[xs,ys,zs] = Geometry3D.FickToSphere(deg2rad(az),deg2rad(el));
[x,y,z] = Geometry3D.FickToSphere(deg2rad([0 50]),deg2rad([0 40]));

dazdx = -y;
dazdy = 1 -  ( y.^2 ) ./ (1 + x);
dazdz = -( y .* z ) ./ (1 + x) ;

deldx = -z ;
deldy = -( y .* z ) ./ (1 + x) ;
deldz = 1 -  ( z.^2 ) ./ (1 + x) ;


rdazwx = 0*dazdx + z.*dazdy - y.*dazdz;
rdazwy = -z.*dazdx - 0.*dazdy + x.*dazdz;
rdazwz = y.*dazdx - x.*dazdy + 0*dazdz;

rdelwx = 0*deldx + z.*deldy - y.*deldz ;
rdelwy = -z.*deldx + 0*deldy + x.*deldz;
rdelwz = y.*deldx - x.*deldy + 0.*deldz ;



figure('color','w')
subplot(1,1,1,'nextplot','add')
axis
% draw sphere
mesh(xs,ys,zs,'FaceAlpha', 1,'facecolor',0.8*[1 1 1]);
xlabel(gca,'x')
ylabel(gca,'y')
zlabel(gca,'z')
% line of sight
line([0 R*1.5 ],[0 0 ],[0 0 ],'color','r','linewidth',2)
line([0 x(2)*1.5 ],[0 y(2)*1.5 ],[0 z(2)*1.5 ],'color','r','linewidth',2)

view(125,15);
axis equal;


hl = [];
hl(1) = quiver3(x,y,z, dazdx, dazdy, dazdz,'linewidth',2,'DisplayName', 'u');
hl(2) = quiver3(x,y,z, deldx, deldy, deldz,'linewidth',2,'DisplayName', 'v');
% quiver3(x,y,z, x, y, z,'linewidth',1)

legend(hl,'Location', 'northeast','box','off','fontsize',14);