%% Visualization of different spherical coordinate systems and motion flows projected onto them
%
% x y z is a right handed reference frame
%
% The eye points towards x and up is z
%
%

clear all, close all

R = 0.025; % radius of the eye
step = 10;
range = 80;
[az, el] = meshgrid(-range:step:range,-range:step:range); % azimuths and elevations to include
coordinateSystems = { 'Listings', 'Fick','Helmholtz', 'Harms','Hess'};


WHICHFLOW = 'linear';
WHICHFLOW = 'rotational';


f = figure('color','white');
f.Position = [f.Position(1) f.Position(2) 4*f.Position(3) f.Position(4)];

for i=1:length(coordinateSystems)
    sys = coordinateSystems{i};
    subplot(4,length(coordinateSystems),i,'nextplot','add');

    % calculate spherical coordinates depending on the coordinate system
    switch(sys)
        case 'Fick'
            [x,y,z] = FickToSphere(az,el);
        case 'Helmholtz'
            [x,y,z] = HelmholtzToSphere(az,el);
        case 'Listings'
            [x,y,z] = ListingsToSphere(az,el);
        case 'Harms'
            [x,y,z] = HarmsToSphere(az,el);
        case 'Hess'
            [x,y,z] = HessToSphere(az,el);
    end

    % scale by radius so the translations are in the right units (m)
    z = z*R;
    x = x*R;
    y = y*R;

    % draw sphere
    view(125,15)
    axis equal;
    mesh(x,y,z,'FaceAlpha', 0.9,'facecolor',[1 1 1]);
    xlabel(gca,'x')
    ylabel(gca,'y')
    zlabel(gca,'z')
    % line of sight
    line([0 R*1.5 ],[0 0 ],[0 0 ],'color','r','linewidth',2)
    title(sys)



    % calculate jacobians for linear or rotational velocity
    % That is, how much a tiny translation or rotation affect the position
    % of a point in the sphere.
    %
    % For rotation is easy because the point remains in the sphere. 
    %
    % For translation is a bit trickier because the point moves outside of
    % the sphere and the vector of motion has to be projected into a
    % tangent plane. 
    %
    % the strategy is to move all the points in the sphere and then
    % calculate the vector from the original point to the new point and
    % project that vector onto the tangent plane by normalizing the second
    % point. So it is a vector between two points in the sphere so we can
    % calcuate it the vector on the actual coordinates of the azimuth and
    % elevation system. For very small displacements the tangent plane and
    % the local surface of the sphere should be approximately the same.

    dv = 0.00001; % differential position change in m
    dw = 0.1; % differential angular rotation in deg

    az2 = az;
    el2 = el;

    colorsflow = {'r' 'b' 'k'};
    dims = {'x' 'y' 'z'};

    dxyz = {[dv,0,0],[0,dv,0],[0,0,dv]}; % assumes distance to the target equal to the radius. It should be scaled by the ratio between the depth and R
    dRM = {RotX(dw), RotY(dw), RotZ(dw)};

    for vi=1:3

        if ( strcmp( WHICHFLOW , 'linear' ) )
            % add a tiny displacement in the corresponding component to all the
            % points in the sphere.
            x2 = x - dxyz{vi}(1);
            y2 = y - dxyz{vi}(2);
            z2 = z - dxyz{vi}(3);
        else
            % add a tiny rotation in the corresponding component to all the
            % points in the sphere.
            pointsrot = [x(:) y(:) z(:)]*dRM{vi};
            x2 = reshape(pointsrot(:,1),size(x));
            y2 = reshape(pointsrot(:,2),size(y));
            z2 = reshape(pointsrot(:,3),size(z));
        end

        % project the tiny displacement onto the sphere
        % (this is the step I am least sure of)
        % But I think it should be correct because it is the projection
        % onto a tangent plane.
        % I do so by simply normalizing the vector so we find the point on
        % the sphere that is closes tot he displaced point. 
        nrm = sqrt(x2.^2 + y2.^2 + z2.^2);
        x2 = x2 ./ nrm;
        y2 = y2 ./ nrm;
        z2 = z2 ./ nrm;

        switch(sys)
            case 'Fick'
                [az2, el2 ] = SphereToFick(x2,y2,z2);
            case 'Helmholtz'
                [az2, el2 ] = SphereToHelmholtz(x2,y2,z2);
            case 'Listings'
                az2 = nan(size(az));
                el2 = nan(size(el));
            case 'Harms'
                [az2, el2 ] = SphereToHarms(x2,y2,z2);
            case 'Hess'
                [az2, el2 ] = SphereToHess(x2,y2,z2);
        end


        subplot(4,length(coordinateSystems),i+length(coordinateSystems)*vi);
        quiver(az, el, az2-az, el2-el, colorsflow{vi});

        switch(WHICHFLOW)
            case 'linear'
                title(['Flow due to ' WHICHFLOW  ' motion along ',dims{vi}])
            case 'rotational'
                title(['Flow due to ' WHICHFLOW  ' motion around ',dims{vi}])
        end
        if ( i > 1)
            xlabel('Azimuth (deg)')
            ylabel('Elevation (deg)')
        else
            xlabel('Angle (deg)')
            ylabel('Eccentricity (deg)')
        end
        grid on

        axis equal;
        set(gca,'xlim',[-80 80],'ylim',[-80 80])
    end

end

% Functions to converte between spherical coordinates and coordinate
% sysstems for 2D rotations

function [x, y, z] = FickToSphere(az, el)
z = sind(el);
y = sind(az).*cosd(el);
x = cosd(el).*cosd(az);
end

function [az,el] = SphereToFick(x,y,z)
el = asind(z);
az = atan2d(y,x);
end

function [x, y, z] = HelmholtzToSphere(az,el)
y = sind(az);
x = cosd(az).*cosd(el);
z = sind(el).*cosd(az);
end

function [az,el] = SphereToHelmholtz(x,y,z)
az = asind(y);
el = atan2d(z,x);
end

function [x, y, z] = ListingsToSphere(angle, eccentricity)
x = cosd(eccentricity);
y = sind(eccentricity).*cosd(angle);
z = sind(eccentricity).*sind(angle);
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

function [az,el] = SphereToHarms(x,y,z)
el = atan2d(z,x);
az = atan2d(y,x);
end

function [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HarmsLinearJacobian(az, el)
    [x, y, z] = HarmsToSphere(az, el);
    dazdx = -y./(y.^2+x.^2);
    dazdy = -y./(y.^2+x.^2);
    dazdz = zeros(size(az));

    deldx = -z./(z.^2+x.^2);
    deldy = zeros(size(el));
    deldz = -z./(z.^2+x.^2);
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

function [az,el] = SphereToHess(~,y,z)
el = asind(z);
az = asind(y);
end



% HORIZONTAL ROTATION right handed
function M = RotZ(theta)
M = [   cos(theta)  -sin(theta)     0;
        sin(theta)  cos(theta)      0;
        0           0               1];
end

% VERTICAL ROTATION right handed
function M = RotY(phi)
M = [   cos(phi)    0               sin(phi);
        0           1               0;
        -sin(phi)   0               cos(phi)];
end

% TORSIONAL ROTATION right handed
function M = RotX(psi)
M = [   1           0               0;
        0           cos(psi)      -sin(psi);
        0           sin(psi)      cos(psi)];
end