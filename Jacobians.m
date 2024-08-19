%%
% close all 
R = 0.025; % radius of the eye
step = 10;
range = 80;
[az, el] = meshgrid(-range:step:range,-range:step:range); % azimuths and elevations to include



[dazdx, dazdy, dazdz, deldx, deldy, deldz] = FickLinearJacobian(az, el);
[dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = FickRotationalJacobian(az, el);

% [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HelmholtzLinearJacobian(az, el);
% [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = HelmholtzRotationalJacobian(az, el);

% [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HarmsLinearJacobian(az, el);
% [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = HarmsRotationalJacobian(az, el);


% [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HessLinearJacobian(az, el);
% [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = HessRotationalJacobian(az, el);

figure('color','w')
colors = get(gca,'colororder');
subplot(2,3,1);
quiver(az, el, dazdx, deldx, 'color', colors(2,:));
axis equal; set(gca,'xlim',[-90 90], 'ylim',[-90 90]); xlabel('azimuzh (deg)'), ylabel('Elevation (deg)'), title('Linear velocity along x = \partial [\theta \psi] / \partial x')
subplot(2,3,2);
quiver(az, el, dazdy, deldy, 'color', colors(1,:));
axis equal; set(gca,'xlim',[-90 90], 'ylim',[-90 90]); xlabel('azimuzh (deg)'), ylabel('Elevation (deg)'), title('Linear velocity along y = \partial [\theta \psi] / \partial y')
subplot(2,3,3);
quiver(az, el, dazdz, deldz, 'color', colors(5,:));
axis equal; set(gca,'xlim',[-90 90], 'ylim',[-90 90]); xlabel('azimuzh (deg)'), ylabel('Elevation (deg)'), title('Linear velocity along z = \partial [\theta \psi] / \partial z')

subplot(2,3,4);
quiver(az, el, dazwxdt, delwxdt, 'color', colors(2,:));
axis equal; set(gca,'xlim',[-90 90], 'ylim',[-90 90]); xlabel('azimuzh (deg)'), ylabel('Elevation (deg)'), title('Angular velocity around x = \partial [\theta \psi] / \omega_xdt')
subplot(2,3,5);
quiver(az, el, dazwydt, delwydt, 'color', colors(1,:));
axis equal; set(gca,'xlim',[-90 90], 'ylim',[-90 90]); xlabel('azimuzh (deg)'), ylabel('Elevation (deg)'), title('Angular velocity around y = \partial [\theta \psi] / \omega_ydt')
subplot(2,3,6);
quiver(az, el, dazwzdt, delwzdt, 'color', colors(5,:));
axis equal; set(gca,'xlim',[-90 90], 'ylim',[-90 90]); xlabel('azimuzh (deg)'), ylabel('Elevation (deg)'), title('Angular velocity around z = \partial [\theta \psi] / \omega_zdt')


figure('color','w')
subplot(1,2,1,'nextplot','add');
colors = get(gca,'colororder');

axis equal; set(gca,'xlim',[-90 90], 'ylim',[-90 90]); xlabel('azimuzh (deg)'), ylabel('Elevation (deg)'), title('Linear Jacobian')
quiver(az, el, dazdx, deldx, 'color', colors(2,:),'LineWidth',1.5, 'DisplayName', '\partial [\theta \psi] / \partial x');
quiver(az, el, dazdy, deldy, 'color', colors(1,:),'LineWidth',1.5, 'DisplayName', '\partial [\theta \psi] / \partial y');
quiver(az, el, dazdz, deldz, 'color', colors(5,:),'LineWidth',1.5, 'DisplayName', '\partial [\theta \psi] / \partial z');

legend('show', 'Location', 'northeast','box','off','fontsize',18);

subplot(1,2,2,'nextplot','add');
axis equal; set(gca,'xlim',[-90 90], 'ylim',[-90 90]); xlabel('azimuzh (deg)'), ylabel('Elevation (deg)'), title('Rotational Jacobian')
quiver(az, el, dazwxdt, delwxdt,'color', colors(2,:),'LineWidth',1.5, 'DisplayName', '\partial [\theta \psi] / \omega_xdt');
quiver(az, el, dazwydt, delwydt,'color', colors(1,:),'LineWidth',1.5, 'DisplayName', '\partial [\theta \psi] / \omega_ydt');
quiver(az, el, dazwzdt, delwzdt,'color', colors(5,:),'LineWidth',1.5, 'DisplayName', '\partial [\theta \psi] / \omega_zdt');

legend('show', 'Location', 'northeast','box','off','fontsize',18);

function [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HarmsLinearJacobian(az, el)
    [x, y, z] = HarmsToSphere(az, el);

    % az = atan2d(y,x);
    % el = atan2d(z,x);
    %
    % x = 1/sqrt(1+tan(az)^2 +tan(el)^2)
    % y = tan(az)*x
    % z = tan(el)*x
    %
    %
    % this are just the derivatives of the SphereToHarms function
    dazdx = rad2deg( -y./(y.^2+x.^2));
    dazdy = rad2deg(x./(y.^2+x.^2));
    dazdz = rad2deg(zeros(size(az)));

    deldx = rad2deg(-z./(z.^2+x.^2));
    deldy = rad2deg(zeros(size(el)));
    deldz = rad2deg(x./(z.^2+x.^2));

    % D = sqrt(1+tand(az).^2+tand(el).^2);
    % 
    % 
    % dazdx = rad2deg( D .* -tand(az)./(1+tand(az).^2));
    % dazdy = rad2deg( D .* 1./(1+tand(az).^2));
    % dazdz = rad2deg(zeros(size(az)));
    % 
    % deldx = rad2deg( D .* -tand(el)./(1+tand(el).^2));
    % deldy = rad2deg(zeros(size(el)));
    % deldz = rad2deg( D .* 1./(1+tand(el).^2));
end

function [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HelmholtzLinearJacobian(az, el)
    [x, y, z] = HarmsToSphere(az, el);
    % this are just the derivatives of the SphereToHarms function
    dazdx = zeros(size(az));
    dazdy = rad2deg( 1./(sqrt(1-y.^2)) ); 
    dazdz = zeros(size(az));

    deldx = rad2deg(-z./(z.^2+x.^2));
    deldy = zeros(size(el));
    deldz = rad2deg(x./(z.^2+x.^2));
end

function [dazdx, dazdy, dazdz, deldx, deldy, deldz] = FickLinearJacobian(az, el)
    [x, y, z] = FickToSphere(az, el);
    % this are just the derivatives of the SphereToHarms function

    % D = sqrt(x.^2+y.^2+z.^2);
    % az = atan2d(y,x);
    % el = asind(z/D);

    dazdx = rad2deg( -y./(y.^2+x.^2));
    dazdy = rad2deg(x./(y.^2+x.^2));
    dazdz = rad2deg(zeros(size(az)));

    DD = sqrt(x.^2+y.^2);
    deldx = -z.*x ./ DD;
    deldy = -z.*y ./ DD;
    deldz =  DD;
end

function [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HessLinearJacobian(az, el)
    [x, y, z] = FickToSphere(az, el);
    % this are just the derivatives of the SphereToHarms function
    dazdx = zeros(size(az));
    dazdy = rad2deg( 1./(sqrt(1-y.^2)) );
    dazdz = zeros(size(az));

    deldx = zeros(size(el));
    deldy = zeros(size(el));
    deldz = rad2deg( 1./(sqrt(1-z.^2)) );
end


function [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = HelmholtzRotationalJacobian(az, el)
    [x, y, z] = HarmsToSphere(az, el);

    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HelmholtzLinearJacobian(az, el);

    % this uses the property that the linear velocity of a point on the
    % sphere rotating by an angular velocity is the cross product of the
    % coordinates of the point and the angular velocity. 
    %
    % so then we can multiply the linear jacobian by the cross product
    % matrix equivalent of the point coordinates

    dazwxdt = 0*dazdx + z.*dazdy - y.*dazdz;
    dazwydt = -z.*dazdx - 0.*dazdy + x.*dazdz;
    dazwzdt = y.*dazdx - x.*dazdy + 0*dazdz;

    delwxdt = 0*deldx + z.*deldy - y.*deldz;
    delwydt = -z.*deldx + 0*deldy + x.*deldz;
    delwzdt = y.*deldx - x.*deldy + 0.*deldz;
end


function [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = HarmsRotationalJacobian(az, el)
    [x, y, z] = HarmsToSphere(az, el);

    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HarmsLinearJacobian(az, el);

    % this uses the property that the linear velocity of a point on the
    % sphere rotating by an angular velocity is the cross product of the
    % coordinates of the point and the angular velocity. 
    %
    % so then we can multiply the linear jacobian by the cross product
    % matrix equivalent of the point coordinates

    dazwxdt = 0*dazdx + z.*dazdy - y.*dazdz;
    dazwydt = -z.*dazdx - 0.*dazdy + x.*dazdz;
    dazwzdt = y.*dazdx - x.*dazdy + 0*dazdz;

    delwxdt = 0*deldx + z.*deldy - y.*deldz;
    delwydt = -z.*deldx + 0*deldy + x.*deldz;
    delwzdt = y.*deldx - x.*deldy + 0.*deldz;
end


function [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = FickRotationalJacobian(az, el)
    [x, y, z] = FickToSphere(az, el);

    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = FickLinearJacobian(az, el);

    % this uses the property that the linear velocity of a point on the
    % sphere rotating by an angular velocity is the cross product of the
    % coordinates of the point and the angular velocity. 
    %
    % so then we can multiply the linear jacobian by the cross product
    % matrix equivalent of the point coordinates

    dazwxdt = 0*dazdx + z.*dazdy - y.*dazdz;
    dazwydt = -z.*dazdx - 0.*dazdy + x.*dazdz;
    dazwzdt = y.*dazdx - x.*dazdy + 0*dazdz;

    delwxdt = 0*deldx + z.*deldy - y.*deldz;
    delwydt = -z.*deldx + 0*deldy + x.*deldz;
    delwzdt = y.*deldx - x.*deldy + 0.*deldz;
end


function [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = HessRotationalJacobian(az, el)
    [x, y, z] = HessToSphere(az, el);

    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HessLinearJacobian(az, el);

    % this uses the property that the linear velocity of a point on the
    % sphere rotating by an angular velocity is the cross product of the
    % coordinates of the point and the angular velocity. 
    %
    % so then we can multiply the linear jacobian by the cross product
    % matrix equivalent of the point coordinates

    dazwxdt = 0*dazdx + z.*dazdy - y.*dazdz;
    dazwydt = -z.*dazdx - 0.*dazdy + x.*dazdz;
    dazwzdt = y.*dazdx - x.*dazdy + 0*dazdz;

    delwxdt = 0*deldx + z.*deldy - y.*deldz;
    delwydt = -z.*deldx + 0*deldy + x.*deldz;
    delwzdt = y.*deldx - x.*deldy + 0.*deldz;
end


function [az,el] = SphereToHarms(x,y,z)
az = atan2d(y,x);
el = atan2d(z,x);
end


function [x, y, z] = HarmsToSphere(az, el)
x = sqrt(1./(1+tand(az).^2+tand(el).^2));
x(az > 90 | az < -90) = -x(az > 90 | az < -90);
x(az == 90 | az == -90 | el == 90 | el == -90) = 0;

z = tand(el).*x;
y = tand(az).*x;

z(el == 90 ) = 1;
z( el == -90) = 1;
y(el == 90 | el == -90) = 0;

z(az == 90 | az == -90) = 0;
y(az == 90 ) = 1;
y(az == -90) = -1;
end


function [x, y, z] = FickToSphere(az, el)
z = sind(el);
y = sind(az).*cosd(el);
x = cosd(el).*cosd(az);
end

function [az,el] = SphereToFick(x,y,z)
D = sqrt(x.^2+y.^2+z.^2);
el = asind(z/D);
az = atan2d(y,x);
end


function [x, y, z] = HelmholtzToSphere(az,el)
y = sind(az);
x = cosd(az).*cosd(el);
z = sind(el).*cosd(az);
end

function [az,el] = SphereToHelmholtz(x,y,z)
D = sqrt(x.^2+y.^2+z.^2);
az = asind(y/D);
el = atan2d(z,x);
end


function [x, y, z] = HessToSphere(az, el)
z = sind(el);
y = sind(az);
x = sqrt(1-z.^2-y.^2);
% Need to flip negative azimuts
x(az >= 90 | az <= -90) = -x(az >= 90 | az <= -90);
% Need to force zero at 90 deg azimuth to avoid some complex
% numbers that can appear numerically
x(az == 90 | az == -90) = 0;

% some of the combinations of azimuth and elevation actually
% don't exist within the sphere.
outsidepoints = (z.^2+y.^2)>=1;
z(outsidepoints) = nan;
y(outsidepoints) = nan;
x(outsidepoints) = nan;
% x = real(x);
end

function [az,el] = SphereToHess(x,y,z)
D = sqrt(x.^2+y.^2+z.^2);
el = asind(z/D);
az = asind(y/D);
end