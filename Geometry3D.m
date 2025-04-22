classdef Geometry3D
    %Geometry3D Helper class for all kinds of 3D operations
    %
    %   Coordinate system
    %
    %       For eye reference frame
    %       Z is positive up
    %       X is positive forward
    %       Y is positive left
    %
    %   Resources
    %
    %       - Haslwanter 1995. Mathematics of three-dimensional eye rotations.
    %       - https://work.thaslwanter.at/thLib/html/rotmat.html

    methods(Static)


        function demoCoordinateSystemsAndPlane()

            app = InteractiveUI('Coordinate systems',@(app) (Geometry3D.demoCoordinateSystemsAndPlaneUpdate(app)), 0.1);
            app.AddDropDown('Coordinate System',   1,  [ "TangentSphere", "Fick", "Helmholtz", "Harms","Hess","ImagePlane"])
            app.AddSlider('Azimuth',           -20,  [-90 90])
            app.AddSlider('Elevation',           20,  [-90 90])
            app.AddSlider('Angular velocity X (deg/s)',   60,  [-100 100] )
            app.AddSlider('Angular velocity Y (deg/s)',   -20,  [-100 100])
            app.AddSlider('Angular velocity Z (deg/s)',   60,  [-100 100])
            app.AddSlider('Linear velocity X (m/s)',    -.5,  [-5 5])
            app.AddSlider('Linear velocity Y (m/s)',    0.8,  [-5 5])
            app.AddSlider('Linear velocity Z (m/s)',    0.7,  [-5 5])
            app.AddDropDown('Stimulus',      2,  ["Ground plane" "Sphere at 1m"])
            app.AddSlider('Height (m)',    1,  [0 10])
            app.AddSlider('Eye elevation (deg)', 0, [-90 90])
            % app.AddDropDown('View3D',             1,  ["Oblique" "TOP" "SIDE"])

            app.AddMenu('Reset to zero', @ResetToZero)

            app.Open();

            function ResetToZero(app)
                app.Values.Azimuth = 0;
                app.Values.Elevation =0;
                app.Values.AngularVelocityX_deg_s_= 0;
                app.Values.AngularVelocityY_deg_s_= 0;
                app.Values.AngularVelocityZ_deg_s_= 0;
                app.Values.LinearVelocityX_m_s_= 0;
                app.Values.LinearVelocityY_m_s_= 0;
                app.Values.LinearVelocityZ_m_s_= 0;
            end
        end

        function [visualDirections, az, el] = SampleVisualDirections(N, type, range)
            % Samples points in the front of a sphere according to the
            % coordinate system according to spiral
            %
            % N: number of visual directions
            % type: type of sample
            %           - 'Fick': Uniform according to Fick coordinates.
            %           - 'Helmholtz': Uniform according to Helmholtz coordinates (.
            %           - 'Harms': Uniform according to Harns coordinates (long, long).
            %           - 'Hess': Uniform according to Hess coordinates (lat, lat).
            %           - 'Spiral': Uniform according to Fick coordinates.
            %           - 'Random': Uniform according to Fick coordinates.
            % range: range of angles to sample from (+- range) in deg


            if ( ~exist('type','var'))
                type = 'Spiral';
            end
            if ( ~exist('range','var'))
                range = 90;
            end

            % get the sample visual directions in spherical coordinates
            % depending on the coordinate system
            N = round(sqrt(N)).^2; % make sure N is a square number
            step = range*2/(sqrt(N)-1);
            [az, el] = meshgrid(deg2rad(-range:step:range),deg2rad(-range:step:range)); % azimuths and elevations to include

            switch(upper(type))
                case 'FICK'
                    [x,y,z] = Geometry3D.FickToSphere(az,el);
                case 'HELMHOLTZ'
                    [x,y,z] = Geometry3D.HelmholtzToSphere(az,el);
                case 'HARMS'
                    [x,y,z] = Geometry3D.HarmsToSphere(az,el);
                case 'HESS'
                    [x,y,z] = Geometry3D.HessToSphere(az,el);
                case 'POLAR'
                    [x,y,z] = Geometry3D.PolarToSphere(az,el);
                case 'SPIRAL'
                    [x, y, z] = Geometry3D.SpiralSphere(N);
                    outofrange = acosd(x) > range;
                    x(outofrange) = [];y(outofrange) = [];z(outofrange) = [];
                    az = []; el = [];
                case {'RANDOM' 'TANGENTSPHERE'}
                    [x, y, z] = Geometry3D.RandomSampleSphere(N);
                    outofrange = acosd(x) > range;
                    x(outofrange) = [];y(outofrange) = [];z(outofrange) = [];
                    az = []; el = [];
            end
            visualDirections = [x(:) y(:) z(:)];
        end

        function [motionField, visualDirections, motionFieldLinear, motionFieldRotational] = CalculateMotionFieldHeadReference(N, v, height, eyeAzimuthHelmholtz, eyeElevationHelmholtz, gain, stim, coordSys)

            % rotation matrix of the eye in head reference
            Reye = Geometry3D.Helm2RotMat(eyeAzimuthHelmholtz,  eyeElevationHelmholtz, 0);
            veye = Reye'*v;

            % Depth along the direction of gaze
            [x,~,z] = Geometry3D.HelmholtzToSphere(eyeAzimuthHelmholtz, eyeElevationHelmholtz);
            d = Geometry3D.CalculateDepthField(eyeAzimuthHelmholtz, eyeElevationHelmholtz);
            D = max( -1/height * (x*sind(eyeElevation) + z*cosd(eyeElevation)) , 0);

            % the eye velocity cancels the linear velocity at the fovea
            % times a gain factor
            [Jw, Jv] = CalculateMotionJacobianFields([0,0,0]);
            w = -gain*Jw'*D*Jv*v;

            % get the retinal motion field
            [motionField, visualDirections, motionFieldLinear, motionFieldRotational] = CalculateMotionField(N, w, veye, [0,0,height], eyeElevationHelmholtz, stim, coordSys);
        end

        function depthsDiopter = CalculateDepthField(visualDirections, eyeOrientationRotMat, eyePositionXYZ, stim)
            % Calculate distances (in diopters) for visual directions

            % Rotate visual directions so they are in world reference
            visualDirections = (eyeOrientationRotMat*visualDirections')';

            eyeheight = eyePositionXYZ(3);

            switch(stim)
                case "Ground plane"
                    depthsDiopter = max( -1/eyeheight * visualDirections(:,3) , 0) ;
                case "Sphere at 1m"
                    depthsDiopter = ones(height(visualDirections),1);
            end
        end

        function [Jw, Jv, J2Dcoord, Jw3Dmotion, Jv3Dmotion] = CalculateMotionJacobianFields(visualDirections, coordSys)
            %   [Jw, Jv, J2Dcoord, Jw3Dmotion, Jv3Dmotion] = ...
            %       CALCULATEMOTIONJACOBIANFIELDS(visualDirections, coordSys)
            %
            %   This function calculates the motion jacobians that connect
            %   the 3D linear and angular velocity with the 2D local motion
            %   at a particular visual direction depending on the
            %   coordinate system used. 
            %
            %   INPUTS:
            %       visualDirections  - A Nx3 matrix specifying visual
            %                           directions in 3D space. 
            %       coordSys          - Coordinate system of the local 2D
            %                           motion ('Fick', 'Helmholtz',
            %                           'Harms', Hess', 'TangentSphere').
            %
            %   OUTPUTS:
            %       Jw              - (2x3xN) Jacobian mapping angular
            %                         motion (e.g., rotational) to local 2D motion.
            %       Jv              - (2x3xN) Jacobian mapping linear
            %                         motion (e.g., translational) to local 2D motion.
            %       J2Dcoord        - (2x3xN) Jacobian for partial derivatives w.r.t. 2D coordinate
            %                         transformations.
            %       Jw3Dmotion      - (2x3xN) Jacobian mapping angular
            %                         motion to local motion in 3D. 
            %       Jv3Dmotion      - (2x3xN) Jacobian mapping linear
            %                         motion to local motion in 3D. 
            %
            %   This is calculated for each of the visual directions given as
            %   an input.
            %   if v is the linear velocity, w is the angular velocity and d
            %   is the depth seen at that visual direction in diopters, then
            %   the local 2D motion is described by the equation
            %
            %       m = d * Jv * v + Jw * w
            %
            %   where the matrices Jv and Jw can be factorized into a matrix
            %   accounting for the 3D motion and another matrix accounting
            %   for the 2D coordinate system projection. 
            %
            %   Jv = J2Dcoord * Jv3Dmotion
            %   Jw = J2Dcoord * Jw3Dmotion
            %

            if ( ~exist('coordSys','var'))
                coordSys = 'TangentSphere';
            end
            
            x = visualDirections(:,1);
            y = visualDirections(:,2);
            z = visualDirections(:,3);

            % 1. Get the jacobian of the 2D coordinate system
            switch(lower(coordSys))
                case 'fick'
                    % daz/dx, daz/dy, daz/dz
                    J2Dcoord(1,1,:) =  -y./(x.^2 + y.^2);
                    J2Dcoord(1,2,:) = x./(x.^2 + y.^2);
                    J2Dcoord(1,3,:) = zeros(size(x));

                    % del/dx, del/dy, del/dz
                    J2Dcoord(2,1,:) = zeros(size(x));
                    J2Dcoord(2,2,:) = zeros(size(x));
                    J2Dcoord(2,3,:) =  1./sqrt(1 - z.^2);

                case 'helmholtz'
                    % daz/dx, daz/dy, daz/dz
                    J2Dcoord(1,1,:) = zeros(size(x));
                    J2Dcoord(1,2,:) =  1./sqrt(1 - y.^2);
                    J2Dcoord(1,3,:) = zeros(size(x));

                    % del/dx, del/dy, del/dz
                    J2Dcoord(2,1,:) =  -z./(z.^2 + x.^2);
                    J2Dcoord(2,2,:) = zeros(size(x));
                    J2Dcoord(2,3,:) = x./(z.^2 + x.^2);

                case 'harms'
                    % daz/dx, daz/dy, daz/dz
                    J2Dcoord(1,1,:)  = -y./(y.^2 + x.^2);
                    J2Dcoord(1,2,:)  = x./(y.^2 + x.^2);
                    J2Dcoord(1,3,:)  = zeros(size(x));

                    % del/dx, del/dy, del/dz
                    J2Dcoord(2,1,:)  = -z./(z.^2 + x.^2);
                    J2Dcoord(2,2,:)  = zeros(size(x));
                    J2Dcoord(2,3,:)  = x./(z.^2 + x.^2);

                case 'hess'
                    denom1 = sqrt( 1 - y.^2 ) ;
                    denom2 = sqrt( 1 - z.^2 ) ;

                    % daz/dx, daz/dy, daz/dz
                    J2Dcoord(1,1,:) = 1./denom1 .* -x.*y;
                    J2Dcoord(1,2,:) = 1./denom1 .* ( x.^2 + z.^2 );
                    J2Dcoord(1,3,:) = 1./denom1 .* -y.*z;

                    % del/dx, del/dy, del/dz
                    J2Dcoord(2,1,:) = 1./denom2 .* -x.*z;
                    J2Dcoord(2,2,:) = 1./denom2 .* -y.*z;
                    J2Dcoord(2,3,:) = 1./denom2 .* ( x.^2 + y.^2 );
                    
                case {'tangentsphere' 'polar'}
                    % daz/dx, daz/dy, daz/dz
                    J2Dcoord(1,1,:) = -y;
                    J2Dcoord(1,2,:) = 1 -  ( y.^2 ) ./ (1 + x);
                    J2Dcoord(1,3,:) = -( y .* z ) ./ (1 + x) ;

                    % del/dx, del/dy, del/dz
                    J2Dcoord(2,1,:) = -z ;
                    J2Dcoord(2,2,:) = -( y .* z ) ./ (1 + x) ;
                    J2Dcoord(2,3,:) = 1 -  ( z.^2 ) ./ (1 + x) ;
            end

            % 2. Get the 3D motion jacobians in 3D space
            Jw3Dmotion = zeros(3,3,length(x));

            Jw3Dmotion(1,1,:) =  0;
            Jw3Dmotion(1,2,:) = -z;
            Jw3Dmotion(1,3,:) =  y; 
            Jw3Dmotion(2,1,:) =  z;
            Jw3Dmotion(2,2,:) =  0;
            Jw3Dmotion(2,3,:) = -x; 
            Jw3Dmotion(3,1,:) = -y;
            Jw3Dmotion(3,2,:) =  x;
            Jw3Dmotion(3,3,:) =  0;

            Jv3Dmotion = zeros(3,3,length(x));

            Jv3Dmotion(1,1,:) = x.^2 - 1;
            Jv3Dmotion(1,2,:) = x.*y;
            Jv3Dmotion(1,3,:) = x.*z; 
            Jv3Dmotion(2,1,:) = x.*y;
            Jv3Dmotion(2,2,:) = y.^2 - 1;
            Jv3Dmotion(2,3,:) = y.*z; 
            Jv3Dmotion(3,1,:) = x.*z;
            Jv3Dmotion(3,2,:) = y.*z;
            Jv3Dmotion(3,3,:) = z.^2 - 1;

            % 3. Get the final motion jacobian in 2D space
            Jv = pagemtimes(J2Dcoord, Jv3Dmotion);
            Jw = pagemtimes(J2Dcoord, Jw3Dmotion);
        end

        function [Jw, Jv, J2Dcoord, Jw3Dmotion, Jv3Dmotion] = CalculateMotionJacobianFieldsNumerical(visualDirections, coordSys)

            if ( ~exist('coordSys','var'))
                coordSys = 'TangentSphere';
            end
            
            J2Dcoord = [];

            N = height(visualDirections);
            
            z = visualDirections(:,3);
            x = visualDirections(:,1);
            y = visualDirections(:,2);


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
            % This this placements actually needs to be multiplied by
            % the ration between the eye radius and the distance of the
            % object seen by this receptor. Right now wea assume ration
            % of 1.
            dw = 0.0001; % differential angular rotation in deg


            dxyz = {[dv,0,0],[0,dv,0],[0,0,dv]}; % differential vectors
            dRM = {Geometry3D.RotX(dw), Geometry3D.RotY(dw), Geometry3D.RotZ(dw)}; % differential rotation matrices
            flows = {'linear', 'rotational'};

                switch(lower(coordSys))
                    case 'fick'
                        [az, el] = Geometry3D.SphereToFick(x,y,z);
                    case 'helmholtz'
                        [az, el] = Geometry3D.SphereToHelmholtz(x,y,z);
                    case {'polar' 'sphere' 'tangentsphere'}
                        [az, el] = Geometry3D.SphereToPolar(x,y,z);
                    case 'harms'
                        [az, el] = Geometry3D.SphereToHarms(x,y,z);
                    case 'hess'
                        [az, el] = Geometry3D.SphereToHess(x,y,z);
                end

            for WHICHFLOW=flows

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
                        pointsrot = [x(:) y(:) z(:)]*dRM{vi}';
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

                    switch(lower(coordSys))
                        case 'fick'
                            [az2, el2 ] = Geometry3D.SphereToFick(x2,y2,z2);
                        case 'helmholtz'
                            [az2, el2 ] = Geometry3D.SphereToHelmholtz(x2,y2,z2);
                        case {'tangentsphere' 'polar'}
                            az2 = nan(size(az));
                            el2 = nan(size(el));
                        case 'harms'
                            [az2, el2 ] = Geometry3D.SphereToHarms(x2,y2,z2);
                        case 'hess'
                            [az2, el2 ] = Geometry3D.SphereToHess(x2,y2,z2);
                    end

                    if ( strcmp( WHICHFLOW , 'linear' ) )

                        Jv(1,vi,:) = (az2-az)/dv;
                        Jv(2,vi,:) = (el2-el)/dv;

                        Jv3Dmotion(1,vi,:) = (x2-x)/dv;
                        Jv3Dmotion(2,vi,:) = (y2-y)/dv;
                        Jv3Dmotion(3,vi,:) = (z2-z)/dv;
                    else
                        Jw3Dmotion(1,vi,:) = (x2-x)/dw;
                        Jw3Dmotion(2,vi,:) = (y2-y)/dw;
                        Jw3Dmotion(3,vi,:) = (z2-z)/dw;

                        Jw(1,vi,:) = (az2-az)/dw;
                        Jw(2,vi,:) = (el2-el)/dw;
                    end
                end
            end
        end

        function [motionField, motionFieldLinear, motionFieldRotational, Jv, Jw, D] = CalculateMotionField(visualDirections, w, v, eyePositionXYZ, eyeOrientationRotMat, stim, coordSys)
            if ( ~exist('coordSys','var'))
                coordSys = 'TangentSphere';
            end

            if ( ~exist('stim','var'))
                stim = 'Ground plane';
            end

            if ( ~exist('eyeOrientationRotMat','var'))
                eyeOrientationRotMat = eye(3);
            end

            [Jw, Jv] = Geometry3D.CalculateMotionJacobianFields(visualDirections, coordSys);

            d = Geometry3D.CalculateDepthField(visualDirections, eyeOrientationRotMat, eyePositionXYZ, stim);
            D = diag(d);

            % for rotational assume that there is nothing when distance is
            % infinite, TODO: make this depend on stimulus configuration
            stimPresent = double(diag(d>0));

            % Calculate the motion field at the visual directions
            motionFieldRotational = stimPresent*squeeze(pagemtimes(Jw, w))';
            motionFieldLinear = D*squeeze(pagemtimes(Jv, v))';
            motionField = motionFieldRotational + motionFieldLinear;
        end

        %%
        function demoCoordinateSystemsAndPlaneUpdate(app)


            w = deg2rad( [app.Values.AngularVelocityX_deg_s_ , app.Values.AngularVelocityY_deg_s_ , app.Values.AngularVelocityZ_deg_s_] )';
            v = [app.Values.LinearVelocityX_m_s_, app.Values.LinearVelocityY_m_s_, app.Values.LinearVelocityZ_m_s_]';
            eyeel = app.Values.EyeElevation_deg_;
            h =  app.Values.Height_m_;
            N = 250;
            N = round(sqrt(N)).^2;
            stim = app.Values.Stimulus;

            coordSys = app.Values.CoordinateSystem;
            if ( isfield(app.Data, 'visualDirections'))
                visualDirections = app.Data.visualDirections;
            else
                visualDirections = Geometry3D.SampleVisualDirections(N,coordSys);
                app.Data.visualDirections = visualDirections;
            end
            [motionField, motionFieldLinear, motionFieldRotational] = Geometry3D.CalculateMotionField(visualDirections, w, v, [0,0,h], eyeel, stim, coordSys);
            x = visualDirections(:,1);
            y = visualDirections(:,2);
            z = visualDirections(:,3);

            if ( ~isfield(app.Data,'hs'))
                % Initialize graphics

                app.Data.hs = struct();

                app.Data.hs.f = figure('color','white');

                app.Data.hs.ax1 = subplot(1,2,1,'nextplot','add');
                colors = get(gca,'colororder');

                % draw sphere
                view(125,15)
                axis equal;
                app.Data.hs.meshSphere = mesh(zeros(2,2),zeros(2,2),zeros(2,2),'FaceAlpha', 0.9,'facecolor',[1 1 1]);
                xlabel(gca,'x')
                ylabel(gca,'y')
                zlabel(gca,'z')
                % line of sight
                set(gca,'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1],'ClippingStyle','rectangle')
                app.Data.hs.TextSys = text(0,0,2, 'hi','FontWeight','bold','HorizontalAlignment','center');

                set(gca,'xtick',[],'ytick',[],'ztick',[])
                set(gca,'visible','off')


                % draw axis
                quiver3(0,0,0,1*2.5, 0, 0,'color',colors(2,:),'linewidth',2);
                quiver3(0,0,0,0, 1*1.5, 0, 0,'color',colors(1,:),'linewidth',2);
                quiver3(0,0,0,0, 0, 1*1.5 ,'color',colors(5,:),'linewidth',2);
                text(2.6,0,0, 'x','fontsize',14,'FontWeight','normal','HorizontalAlignment','center');
                text(0,1.7,0, 'y','fontsize',14,'FontWeight','normal','HorizontalAlignment','center');
                text(0,0,1.55, 'z','fontsize',14,'FontWeight','normal','HorizontalAlignment','center');


                % draw velocity vectors
                app.Data.hs.quiverw = quiver3(0,0,0, 1 , -0.5 , 1.3 ,'color',colors(3,:),'linewidth',2,'LineStyle','-.');
                app.Data.hs.textw = text(1,-0.5,1.3, '\omega','fontsize',14,'FontWeight','normal','HorizontalAlignment','right');
                app.Data.hs.quiverv = quiver3(0,0,0, 0.5, 1,-1 ,'color',colors(4,:),'linewidth',2,'LineStyle','-.');
                app.Data.hs.textv = text(0.5,1.1,-1, 'v','fontsize',14,'FontWeight','normal','HorizontalAlignment','left');

                % draw point
                app.Data.hs.meshpoint = mesh(zeros(2,2), zeros(2,2),zeros(2,2), 'EdgeColor','none','FaceColor','k');
                app.Data.hs.textpoint = text(0*1.2,0*1.2,0*1.2, 'p','fontsize',14, 'FontWeight','normal','HorizontalAlignment','right','VerticalAlignment','top');

                % draw tangent plane in 3D
                app.Data.hs.meshtangent = mesh(zeros(2,2), zeros(2,2),zeros(2,2), 'EdgeColor',[0.7 0.7 0.7],'FaceAlpha', 0.3,'facecolor',[0.8 0.8 0.8]);
                app.Data.hs.quiverdaz = quiver3(0,0,0, 0*2, 0*2 ,0*2,'color',[0.6 0.6 0.6],'linewidth',2);
                app.Data.hs.quiverdel = quiver3(0,0,0, 0*2, 0*2 ,0*2,'color',[0.6 0.6 0.6],'linewidth',2);
                app.Data.hs.quivertvw = quiver3(0,0,0, 0*2, 0*2 ,0*2,'color',colors(3,:),'linewidth',2  ,'LineStyle',':');
                app.Data.hs.quivertvv = quiver3(0,0,0, 0*2, 0*2 ,0*2,'color',colors(4,:),'linewidth',2  ,'LineStyle',':');
                app.Data.hs.quivertv  = quiver3(0,0,0, 0*2, 0*2 ,0*2,'color','k','linewidth',2 );


                % draw flat plane
                app.Data.hs.ax2 = subplot(1,2,2,'nextplot','add');
                axis equal;
                % xlabel('Azimuth ')
                % ylabel('Elevantion ')

                set(gca,'xlim',[-90 90],'ylim',[-80 80])
                set(gca,'xticklabels',[])
                set(gca,'yticklabels',[])


                app.Data.hs.textpointflat = text(0*1.2,0*1.2, 'p','fontsize',14, 'FontWeight','normal','HorizontalAlignment','right','VerticalAlignment','top');
                app.Data.hs.quiverdazflat = quiver(0,0, 0*2 ,0*2,'color',[0.6 0.6 0.6],'linewidth',2);
                app.Data.hs.quiverdelflat = quiver(0,0, 0*2 ,0*2,'color',[0.6 0.6 0.6],'linewidth',2 );


                app.Data.hs.quivertJwflat = quiver(0,0, 0*2 ,0*2,'color',colors(3,:),'linewidth',1,'LineStyle',':');
                app.Data.hs.quivertJvflat = quiver(0,0, 0*2 ,0*2,'color',colors(4,:),'linewidth',1,'LineStyle',':');
                app.Data.hs.quivertAllvflat  = quiver(0,0, 0*2 ,0*2,'color','k','linewidth',1);

                legend([app.Data.hs.quivertJvflat app.Data.hs.quivertJwflat app.Data.hs.quivertAllvflat app.Data.hs.quiverdazflat],{'Linear motion' 'Rotational motion' 'Total motion' 'Tanget basis'})

            end

            set(app.Data.hs.TextSys, 'String',coordSys)

            % get data and update


            % get points to draw meshgrid, different from the sampling of
            % the motion field on the sphere

            range = 80;
            step = range*2/(sqrt(N)-1);
            [az, el] = meshgrid(deg2rad(-range:step:range),deg2rad(-range:step:range)); % azimuths and elevations to include
            switch(coordSys)
                case 'Fick'
                    [xs,ys,zs] = Geometry3D.FickToSphere(az,el);
                case 'Helmholtz'
                    [xs,ys,zs] = Geometry3D.HelmholtzToSphere(az,el);
                case 'Harms'
                    [xs,ys,zs] = Geometry3D.HarmsToSphere(az,el);
                case 'Hess'
                    [xs,ys,zs] = Geometry3D.HessToSphere(az,el);
                case 'TangentSphere'
                    [xs,ys,zs] = Geometry3D.FickToSphere(az,el);
                    % change az and el for the flat plot for now in this
                    % simple way
                    az = y.*acos(x);
                    el = z.*acos(x);
            end


            % example point
            paz = deg2rad(app.Values.Azimuth);
            pel = deg2rad(app.Values.Elevation);


            [px,py,pz] = Geometry3D.AzElToSphere(paz,pel, coordSys);
            [Jwpoint, Jvpoint] = Geometry3D.CalculateMotionJacobianFields([px,py,pz],coordSys);

            % switch(coordSys)
            %     case 'Fick'
            %         [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.FickInverseJacobian(paz, pel);
            %         [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.FickCoordinatesJacobian(px,py,pz);
            % 
            %     case 'Helmholtz'
            %         [px,py,pz] = Geometry3D.HelmholtzToSphere(paz,pel);
            %         [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.HelmholtzInverseJacobian(paz, pel);
            %         [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HelmholtzCoordinatesJacobian(px,py,pz);
            %     case 'Harms'
            %         [px,py,pz] = Geometry3D.HarmsToSphere(paz,pel);
            %         [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.HarmsInverseJacobian(paz, pel);
            %         [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HarmsCoordinatesJacobian(px,py,pz);
            %     case 'Hess'
            %         [px,py,pz] = Geometry3D.HessToSphere(paz,pel);
            %         [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.HessInverseJacobian(paz, pel);
            %         [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HessCoordinatesJacobian(px,py,pz);
            %     case 'TangentSphere'
            %         [px,py,pz] = Geometry3D.FickToSphere(paz,pel);
            %         [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.TangentSphereJacobian(px,py,pz);
            %         dxdaz = dazdx;
            %         dydaz = dazdy;
            %         dzdaz = dazdz;
            %         dxdel = deldx;
            %         dydel = deldy;
            %         dzdel = deldz;
            %         [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.TangentRotationalJacobian(px,py,pz);
            % end

            % udpate sphere
            set(app.Data.hs.meshSphere, 'xdata',xs,'ydata',ys,'zdata',zs)


            set(app.Data.hs.quiverw, 'UData',w(1),'VData',w(2),'WData',w(3));
            set(app.Data.hs.quiverv, 'UData',v(1),'VData',v(2),'WData',v(3));
            set(app.Data.hs.textw,'Position',w);
            set(app.Data.hs.textv,'Position',v);

            % update point
            [psx,psy,psz] = sphere(10);
            set(app.Data.hs.meshpoint, 'xdata',psx*0.05+px,'ydata',psy*0.05+py,'zdata',psz*0.05+pz)
            set(app.Data.hs.textpoint, 'Position', [px*1.2,py*1.2,pz*1.2])

            % draw tangent plane
            set(app.Data.hs.quiverdaz, 'xdata',px,'ydata',py,'zdata',pz);
            set(app.Data.hs.quiverdel, 'xdata',px,'ydata',py,'zdata',pz);
            set(app.Data.hs.quiverdaz, 'UData',dxdaz,'VData',dydaz,'WData',dzdaz);
            set(app.Data.hs.quiverdel, 'UData',dxdel,'VData',dydel,'WData',dzdel);

            [daz,del] = meshgrid(-0.5:0.1:1.2,-0.5:0.1:1.2);
            set(app.Data.hs.meshtangent, 'xdata',daz*dxdaz + del*dxdel + px,'ydata',daz*dydaz + del*dydel +py,'zdata',daz*dzdaz + del*dzdel +pz)

            Ji =[dxdaz, dydaz, dzdaz; dxdel, dydel, dzdel]'; % base of the tangent plane
            tvw = Ji*[rdazwx, rdazwy, rdazwz; rdelwx, rdelwy, rdelwz]*w;
            tvv = Ji*[dazdx, dazdy, dazdz; deldx, deldy, deldz]*v;
            tv = tvw + tvv;

            set(app.Data.hs.quivertvw, 'xdata',px,'ydata',py,'zdata',pz);
            set(app.Data.hs.quivertvv, 'xdata',px,'ydata',py,'zdata',pz);
            set(app.Data.hs.quivertv, 'xdata',px,'ydata',py,'zdata',pz);
            set(app.Data.hs.quivertvw, 'UData',tvw(1),'VData',tvw(2),'WData',tvw(3));
            set(app.Data.hs.quivertvv, 'UData',tvv(1),'VData',tvv(2),'WData',tvv(3));
            set(app.Data.hs.quivertv, 'UData',tv(1),'VData',tv(2),'WData',tv(3));



            % update flat motion field

            azdeg = rad2deg(az);
            eldeg = rad2deg(el);
            motionFieldRotational = rad2deg(motionFieldRotational);
            motionFieldLinear = rad2deg(motionFieldLinear);
            motionField = rad2deg(motionField);
            motionFieldRotational(abs(azdeg)>80 | abs(eldeg)>80) = nan;
            motionFieldLinear(abs(azdeg)>80 | abs(eldeg)>80) = nan;
            motionField(abs(azdeg)>80 | abs(eldeg)>80) = nan;

            set(app.Data.hs.quivertJwflat, 'xdata', azdeg(:), 'ydata', eldeg(:), 'UData', motionFieldRotational(:,1), 'VData', motionFieldRotational(:,2));
            set(app.Data.hs.quivertJvflat, 'xdata', azdeg(:), 'ydata', eldeg(:), 'UData', motionFieldLinear(:,1), 'VData', motionFieldLinear(:,2));
            set(app.Data.hs.quivertAllvflat, 'xdata', azdeg(:), 'ydata', eldeg(:), 'UData', motionField(:,1), 'VData', motionField(:,2));



            % update flat motion field point

            pointmotionFieldRotational = [rdazwx, rdazwy, rdazwz; rdelwx, rdelwy, rdelwz] * w;
            pointmotionFieldLinear = [dazdx, dazdy, dazdz; deldx, deldy, deldz] * v;
            pointmotionField = pointmotionFieldRotational + pointmotionFieldLinear;


            pazdeg = rad2deg(paz);
            peldeg = rad2deg(pel);

            pointmotionFieldRotational = rad2deg(pointmotionFieldRotational);
            pointmotionFieldLinear = rad2deg(pointmotionFieldLinear);
            pointmotionField = rad2deg(pointmotionField);

            set(app.Data.hs.textpointflat, 'Position', [pazdeg, peldeg]);
            set(app.Data.hs.quiverdazflat, 'xdata', pazdeg, 'ydata', peldeg);
            set(app.Data.hs.quiverdelflat, 'xdata', pazdeg, 'ydata', peldeg);
            set(app.Data.hs.quiverdazflat, 'UData', rad2deg(1), 'VData', rad2deg(0));
            set(app.Data.hs.quiverdelflat, 'UData', rad2deg(0), 'VData', rad2deg(1));


            drawnow limitrate
        end

    end

    methods(Static) % DEMO DISPARITY

        function demoDisparity()

            helptext = {...
                '--------------------------------------------------------'...
                'Coordinate system and conventions:' ...
                '--------------------------------------------------------'...
                'Right hand cartesian coordinate system.' ...
                'X points forward from the head.'...
                'Y points left from the head.'...
                'Z points up from the head'...
                'Right hand rule for rotations.' ...
                'For some angles we use Helmholtz coordinate system.' ...
                'Azimuth is the horizontal rotation (positive to the right)' ...
                'Elevation is the vertical rotation (positive up)' ...
                '--------------------------------------------------------'...
                'Options:' ...
                '--------------------------------------------------------'...
                'IPD mm : distance between the eye centers in mm.'...
                ''...
                'Stimulus Size cm: size of the stimulus in cm.'...
                'Stimulus Distance cm: distance from the eyes to the stimulus.'...
                'Stimulus slant deg: slant of the stimulus surface in deg. Angle between X axis and surface normal. '...
                'Stimulus Tilt deg: tilt of the stimulus surface in deg. Angle between projection of surface normal into ZY plane and vertical axis.'...
                ''...
                'Fixation Distance cm: Distance from the eyes to the fixation spot.'...
                'Fixation Azimuth deg: horizontal angle in Helmholtz coordinates from the eyes to the fixation spot.'...
                'Fixation Elevation deg: vertical angle in Helmholtz coordinates from the eyes to the fixation spot.'...
                ''...
                'Listings Plane Pitch deg: binocular rotation of the listing''s plane along the y axis. Positive means pitch down.'...
                'Listings L2 factor: how many degrees does listing''s plane rotate along the Z axis (yaw) for reach degree of vergence. Positive causes a temporal rotation of the primary position.' ...
                'Torsion Version deg: cycloversion off of Listing''s plane. Positive means top pole to right shoulder.'...
                'Torsion Vergence deg: cyclovergence off of Listing''s plane. Positive means top pole away the nose (extorsion). Each eye rotates half this ammount.'...
                'Retinal Shear deg (not confirmed accuracy): elevation dependent horizontal retinal shear. Two times angle of the retinal vertical meridian relative to the geometrical meridian. Positive means top pole to the nose (extorsion). For each eye the tilt is half of this ammount'...
                'Hering–Hillebrand deviation (not confirmed accuracy): deviation between the horizontal empirical and geometrical horopters.'...
                ''...
                'View3D: different view points for the 3D plots.'...
                };

            updateIntervalMs = 0.2;

            app = InteractiveUI('Disparity Simulator',@(app) (Geometry3D.demoDisparityUpdate(app)), updateIntervalMs, helptext);
            app.AddDropDown('Stimulus',                 1,  ["CROSS" "FIXATION" "RANDOMPLANE" "GRIDPLANE" "HLINE" "VLINE"])
            app.AddSlider('IPD mm',                     60, [10 100])
            app.AddSlider('Stimulus Size cm',           40, [10 200])
            app.AddSlider('Stimulus Distance cm',       40, [10 200])
            app.AddSlider('Stimulus slant deg',         0,  [-90 90])
            app.AddSlider('Stimulus Tilt deg',          0,  [0 90])
            app.AddSlider('Fixation Distance cm',       30, [10 200])
            app.AddSlider('Fixation Azimuth deg',       0,  [-50 50])
            app.AddSlider('Fixation Elevation deg',     0,  [-50 50])
            app.AddSlider('Listings Plane Pitch deg',   0,  [-20 20])
            app.AddSlider('Listings L2 factor',         0,  [-2 2])
            app.AddSlider('Torsion Version deg',        0,  [-20 20])
            app.AddSlider('Torsion Vergence deg',       0,  [-10 10])
            app.AddSlider('Retinal Shear deg',          0,  [-5 5])
            app.AddSlider('Hering–Hillebrand deviation',0,  [-0.5 0.5])
            app.AddDropDown('View3D',                   1,  ["OBLIQUE" "TOP" "BACK" "SIDE"])
            app.AddSlider('Screen slant deg',           0,  [-30 30])
            app.AddDropDown('Extended horopter',                   1,  ["YES" "NO"])


            app.Data.Screen = struct();
            app.Data.Screen.SizeCm = [30*16/9 30];
            app.Data.Screen.ResPix = [1920 1080];
            app.Data.Screen.DistanceCm = 57;
            app.Data.Screen.SlantDeg = 0;

            app.Data.FixationSpot = struct();
            app.Data.FixationSpot = [0 0 0]';

            app.Open();

        end

        function h = demoDisparityInitPlots()

            % struct containing all the graphical objects
            h = struct();

            % create the figure
            screen_size = get(0,'ScreenSize');
            margin = floor(0.1*(screen_size(4)));
            h.figure = figure('color','w','position',floor([margin margin screen_size(3)*2.8/4 screen_size(4)*2/4 ]));

            plot3Dposition              = [-0.03       0    0.7    1.5];
            plotRetinaPosition          = [0.44    0.3    0.28    0.7];
            plotDisparityPosition       = [0.71    0.3    0.3    0.7];
            plotHaploscopeLeftPosition  = [0.65    0.01    0.15    0.35];
            plotHaploscopeRightPosition = [0.80    0.01    0.15    0.35];
            plotListingsPosition        = [0.5     0.01    0.15    0.35];

            %-----------------------------
            % 3D plot
            %-----------------------------
            h.plot3D.ax = axes('OuterPosition', plot3Dposition, 'nextplot','add');


            grid
            ylim([-100 100]), zlim([-100 100]), xlim([-5 200]);
            xlabel('X (cm)'), zlabel('Z (cm)'), ylabel('Y (cm)');
            % view(-60,10)
            title('3D world');

            % world points, shadow, and fixation
            h.plot3D.wordlPoints        = line(0,0, 0,'linestyle','none','marker','o','Color','k');
            h.plot3D.wordlPointsShadow  = line(0,0, 0,'linestyle','none','marker','o','Color',[0.8 0.8 0.8]);
            h.plot3D.fixationSpot       = line(0,0, 0,'linestyle','none','marker','o','Color','r','LineWidth',2, 'markersize',20);

            % eyes
            h.plot3D.LeftEye.Center     = plot3(0, 0, 0,'o','Color','b', 'markersize',10); % left eye fixation spot and right eye
            h.plot3D.RightEye.Center    = plot3(0, 0, 0,'o','Color','r', 'markersize',10); % left eye fixation spot and right eye
            h.plot3D.LeftEye.Lines      = line(0,0,0,'linestyle','-','Color','b'); 
            h.plot3D.RightEye.Lines     = line(0,0,0,'linestyle','-','Color','r'); 

            % horopter surface
            h.plot3D.horopter = surf(zeros(2),zeros(2),zeros(2),zeros(2));
            originalCmap = colormap('turbo');
            set(h.plot3D.horopter,'EdgeColor','none');
            set(h.plot3D.horopter,'FaceAlpha',0.7);
            set(h.plot3D.ax,'Projection','perspective')
            set(h.plot3D.ax,'CameraPosition',[-100 -100 100],'cameratarget',[200 0 0])

            nColors = size(originalCmap, 1);
            % Create a normalized index vector from 0 to 1
            xi = linspace(0, 1, nColors);
            xi2 = linspace(0, 1, nColors*4);

            % Apply a transformation to compact one side.
            % Here, squaring the indices (xi.^2) compresses values near 0 more than near 1.
            xi_transformed = xi2.^0.5;  % Adjust the exponent as needed

            % Interpolate to create the new colormap based on the transformed indices
            newCmap = interp1(xi, originalCmap, xi_transformed);
            colormap(newCmap);
            set(gca,'clim',[0 50])

            % legend for 3d plot
            h.plot3D.legend = legend([h.plot3D.LeftEye.Lines, h.plot3D.RightEye.Lines, h.plot3D.fixationSpot, h.plot3D.wordlPoints], {'Left eye','Right eye','Fixation point','Stimulus'},'Location','southwest','box','off');

            %-----------------------------
            % retina
            %-----------------------------
            
            h.plotRetina.ax = axes('OuterPosition', plotRetinaPosition, 'nextplot','add');
            title({'Retinal image - visual directions' '(equi-eccentricity)'});

            axis equal; % Ensure equal scaling on both axes
            set(gca,"XLim",[-90 90],"YLim",[-90 90])
            grid off; axis off;
            set(gca,'clim',[0 50])
            c = colorbar;
            c.Label.String = 'Min Disparity (deg)';

            % Define parameters for circles and radial lines
            R = 90;
            numCircles = 6;                   % Number of concentric circles
            numRadialLines = 12;              % Number of radial lines
            labelAngle = pi/8;                % Position labels at 45 degrees for clarity
            Geometry3D.DrawPolarGrid(R, numCircles, numRadialLines, labelAngle);

            h.plotRetina.horopter = surf(zeros(2),zeros(2),zeros(2),zeros(2));
            colormap(newCmap);
            h.plotRetina.horopter.EdgeColor = 'none';
            h.plotRetina.horopter.FaceAlpha = 0.4;

            % stimulus in the left and right eye
            h.plotRetina.leftEyePoints= plot(0, 0,'bo');
            h.plotRetina.rightEyePoints = plot(0, 0,'ro');

            % fixation point
            h.plotRetina.fixationSpot = plot(0,0,'ro','linewidth',3); % TODO: add left and right eye fixation point for missaligned eyes
            legend([h.plotRetina.leftEyePoints, h.plotRetina.rightEyePoints , h.plotRetina.fixationSpot], {'Left eye','Right eye','Fixation point'},'Location','northeast','box','off')


            %-----------------------------
            % disparity plot
            %-----------------------------

            h.plotDisparity.ax = axes('OuterPosition', plotDisparityPosition, 'nextplot','add');
            set(gca, 'PlotBoxAspectRatio',[1 1 1])
            grid
            set(gca,'xlim',[-40 40],'ylim',[-40 40])
            xlabel('Horizontal disparity (Helmtholz, deg, positive uncrossed)')
            ylabel('Vertical disparity (Helmtholz, deg)')

            h.plotDisparity.disparities = quiver(0,0, 0, 0, 'AutoScale', "off", 'MaxHeadSize', 0.1);

            %-----------------------------
            % haploscope plot
            %-----------------------------

            h.plotHaploscope = [];

            h.plotHaploscope.axLeft = axes('OuterPosition', plotHaploscopeLeftPosition, 'nextplot','add');
            h.plotHaploscope.leftPoints     = plot(0, 0, 'bo');
            h.plotHaploscope.leftFixation   = plot(0, 0, 'ro', 'linewidth',2,'markersize',15);
            grid
            set(gca,'xlim',[0 1920],'ylim',[0 1080])
            set(gca,'PlotBoxAspectRatio',[16 9 1])
            xlabel({'Haploscope render left eye screen' '53.3cm x 30 cm, 1920x80, 57 cm dist'});
            set(gca,'box','on')

            h.plotHaploscope.axRight = axes('OuterPosition',plotHaploscopeRightPosition, 'nextplot','add');
            h.plotHaploscope.rightPoints  = plot(0, 0, 'ro');
            h.plotHaploscope.rightFixation  = plot(0, 0, 'ro', 'linewidth',2,'markersize',15);
            grid
            set(gca,'xlim',[0 1920],'ylim',[0 1080])
            set(gca,'yticklabels',[])
            set(gca,'PlotBoxAspectRatio',[16 9 1])
            xlabel({'Haploscope render right eye screen' '53.3cm x 30 cm, 1920x80, 57 cm dist'});
            set(gca,'box','on')

            %-----------------------------
            % Listing planes plot
            %-----------------------------

            % eyes
            h.plotListings.ax = axes('OuterPosition',plotListingsPosition, 'nextplot','add');
            view(h.plotListings.ax, -40 ,  40)
            axis equal
            grid
            set(h.plotListings.ax,'xlim',[-15 15],'ylim',[-15 15],'zlim',[-15 15])
            title('Listing''s planes')

            % eyes
            h.plotListings.LeftEye.Center     = plot3(0, 0, 0,'o','Color','b', 'markersize',10); % left eye fixation spot and right eye
            h.plotListings.RightEye.Center    = plot3(0, 0, 0,'o','Color','r', 'markersize',10); % left eye fixation spot and right eye
            h.plotListings.LeftEye.Lines      = line(0,0,0,'linestyle','-','Color','b'); 
            h.plotListings.RightEye.Lines     = line(0,0,0,'linestyle','-','Color','r'); 

            % Listing planes surface
            h.plotListings.LeftPlane    = surf(zeros(2),zeros(2),zeros(2),zeros(2),'facecolor','white','edgecolor','b');
            h.plotListings.RightPlane   = surf(zeros(2),zeros(2),zeros(2),zeros(2),'facecolor','white','edgecolor','r');
        end

        function demoDisparityUpdate(app)

            Values = app.Values;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % make the world
            worldPoints = Geometry3D.MakeWorldPoints(Values.Stimulus, Values.StimulusSizeCm, Values.StimulusDistanceCm, Values.StimulusTiltDeg, Values.StimulusSlantDeg);

            % Make the eyes
            eyes = Geometry3D.MakeEyes(Values.IPDMm/10, app.Data.FixationSpot, Values.TorsionVersionDeg, Values.TorsionVergenceDeg, Values.ListingsPlanePitchDeg, Values.ListingsL2Factor, Values.RetinalShearDeg, Values.Hering_HillebrandDeviation);
            horopter = Geometry3D.ComputeExtendedHoropter(eyes);

            % update the fixation spot. Note the flip of x and add it to
            % the world points to convert to eye and screen points.
            [fx, fy, fz] = Geometry3D.HelmholtzToSphere( deg2rad(-Values.FixationAzimuthDeg), deg2rad(Values.FixationElevationDeg) );
            app.Data.FixationSpot = Values.FixationDistanceCm * [fx, fy, fz]';
            worldPoints{end+1,{'X' 'Y' 'Z'}} = app.Data.FixationSpot';

            % setup the screens of the haploscope
            leftEyeScreen = Geometry3D.MakeScreen(eyes.L.Center, app.Data.Screen.DistanceCm, app.Data.Screen.SizeCm, app.Data.Screen.ResPix, Values.ScreenSlantDeg);
            rightEyeScreen = Geometry3D.MakeScreen(eyes.R.Center,app.Data.Screen.DistanceCm, app.Data.Screen.SizeCm, app.Data.Screen.ResPix, Values.ScreenSlantDeg);
           
            
            %% Get the points in eye and screen coordinates
            eyePoints = Geometry3D.Points3DToEyes(worldPoints, eyes);
            leftScreenPoints = Geometry3D.PointsEyesToScreen( worldPoints, leftEyeScreen);
            rightScreenPoints = Geometry3D.PointsEyesToScreen( worldPoints, rightEyeScreen);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update graphics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if ( ~isfield(app.Data, "h") || ~isfield(app.Data.h, "figure") || ~isvalid(app.Data.h.figure))
                % If figure does not exist create it with all the plots and
                % handles to them
                app.Data.h = Geometry3D.demoDisparityInitPlots();
            end

            h = app.Data.h;

            % -----------------------------------
            % update 3D plot
            % -----------------------------------

            set(h.plot3D.fixationSpot,      'xdata', app.Data.FixationSpot(1), 'ydata', app.Data.FixationSpot(2), 'zdata', app.Data.FixationSpot(3));

            set(h.plot3D.LeftEye.Lines,     'xdata', eyes.L.LinesForDrawing(:,1),   'ydata', eyes.L.LinesForDrawing(:,2), 'zdata', eyes.L.LinesForDrawing(:,3));
            set(h.plot3D.LeftEye.Center,    'xdata', eyes.L.Center(1),              'ydata', eyes.L.Center(2), 'zdata', eyes.L.Center(3));
            set(h.plot3D.RightEye.Lines,    'xdata', eyes.R.LinesForDrawing(:,1),   'ydata', eyes.R.LinesForDrawing(:,2), 'zdata', eyes.R.LinesForDrawing(:,3));
            set(h.plot3D.RightEye.Center,   'xdata', eyes.R.Center(1),              'ydata', eyes.R.Center(2), 'zdata', eyes.R.Center(3));

            set(h.plot3D.wordlPoints,       'xdata', worldPoints.X, 'ydata',worldPoints.Y, 'zdata', worldPoints.Z);
            set(h.plot3D.wordlPointsShadow, 'xdata', worldPoints.X, 'ydata',worldPoints.Y, 'zdata', h.plot3D.ax.ZLim(1)*ones(size(worldPoints.Z)));

            h3d = horopter.Points;
            h3d(abs(horopter.Az) > deg2rad(60) | abs(horopter.El) > deg2rad(60)) = nan; % restrict the eccentricities
            if ( app.Values.ExtendedHoropter == "YES")
                set(h.plot3D.horopter, 'xdata', h3d(:,:,1), 'ydata', h3d(:,:,2), 'zdata', h3d(:,:,3), 'cdata', horopter.MinDisparity);
            else
                set(h.plot3D.horopter, 'xdata', 0, 'ydata', 0, 'zdata', 0, 'cdata', 0);
            end

            if (~isfield(app.Data, 'View3D') || app.Data.View3D ~= string(app.Values.View3D) )
                switch(app.Values.View3D)
                    case 'OBLIQUE'
                        set(h.plot3D.ax,'OuterPosition', [-0.03       0    0.7    1.5]);
                        set(h.plot3D.ax,'Projection','perspective')
                        set(h.plot3D.ax,'CameraPosition',[-100 -100 100],'cameratarget',[200 0 0])
                    case 'TOP'
                        set(h.plot3D.ax,'OuterPosition', [0       0    0.5    1]);
                        set(h.plot3D.ax,'Projection','orthographic')
                        view([h.plot3D.ax h.plotListings.ax], 0, 90)
                    case 'BACK'
                        set(h.plot3D.ax,'OuterPosition', [0       0    0.5    1]);
                        set(h.plot3D.ax,'Projection','orthographic')
                        view([h.plot3D.ax h.plotListings.ax], -90, 0)
                    case 'SIDE'
                        set(h.plot3D.ax,'OuterPosition', [0       0    0.5    1]);
                        set(h.plot3D.ax,'Projection','orthographic')
                        view([h.plot3D.ax h.plotListings.ax], 0, 0)
                end
                app.Data.View3D = string(Values.View3D);
            end

            % -----------------------------------
            % update retina plot
            % -----------------------------------
            
            % get azimuthal equidistant projections for horopter surface
            % and stimulus points
            ph = rad2deg(Geometry3D.AzimuthalEquidistantProjectionZY(horopter.VisualDirections));
            l = rad2deg(Geometry3D.AzimuthalEquidistantProjectionZY(eyePoints{:,{'LX' 'LY' 'LZ'  }}));
            r = rad2deg(Geometry3D.AzimuthalEquidistantProjectionZY(eyePoints{:,{'RX' 'RY' 'RZ'  }}));

            ph = reshape(ph,size(horopter.Points));
            if ( app.Values.ExtendedHoropter == "YES")
                set(h.plotRetina.horopter,          'xdata', -ph(:,:,2), 'ydata', ph(:,:,3), 'zdata', 0*ph(:,:,2), 'cdata', horopter.MinDisparity);
            else
                set(h.plotRetina.horopter,          'xdata', 0, 'ydata',0, 'zdata', 0, 'cdata', 0);
            end
            set(h.plotRetina.leftEyePoints,     'XData', -l(:,2), 'YData',l(:,3));
            set(h.plotRetina.rightEyePoints,    'XData', -r(:,2), 'YData',r(:,3));

            % update disparity plot
            set(h.plotDisparity.disparities, ...
                'xdata',-eyePoints.RH, 'ydata', eyePoints.RV, ...
                'udata', eyePoints.HDisparity, 'vdata', eyePoints.VDisparity);

            % update screen plots
            set(h.plotHaploscope.leftPoints, 'xdata', leftScreenPoints(1:end-1, 1), 'ydata',leftScreenPoints(1:end-1, 2));
            set(h.plotHaploscope.rightPoints, 'xdata', rightScreenPoints(1:end-1, 1), 'ydata',rightScreenPoints(1:end-1, 2));
            set(h.plotHaploscope.leftFixation, 'xdata', leftScreenPoints(end, 1), 'ydata', leftScreenPoints(end, 2));
            set(h.plotHaploscope.rightFixation, 'xdata', rightScreenPoints(end, 1), 'ydata', rightScreenPoints(end, 2));

            % update listing's plot
            set(h.plotListings.LeftEye.Lines, 'xdata', eyes.L.LinesForDrawing(:,1), 'ydata', eyes.L.LinesForDrawing(:,2), 'zdata', eyes.L.LinesForDrawing(:,3));
            set(h.plotListings.RightEye.Lines, 'xdata', eyes.R.LinesForDrawing(:,1), 'ydata', eyes.R.LinesForDrawing(:,2), 'zdata', eyes.R.LinesForDrawing(:,3));
            set(h.plotListings.LeftEye.Center, 'xdata', eyes.L.Center(1), 'ydata', eyes.L.Center(2), 'zdata', eyes.L.Center(3));
            set(h.plotListings.RightEye.Center, 'xdata', eyes.R.Center(1), 'ydata', eyes.R.Center(2), 'zdata', eyes.R.Center(3));
            set(h.plotListings.LeftPlane, 'xdata', eyes.L.ListingPlaneX, 'ydata', eyes.L.ListingPlaneY, 'zdata', eyes.L.ListingPlaneZ);
            set(h.plotListings.RightPlane, 'xdata', eyes.R.ListingPlaneX, 'ydata', eyes.R.ListingPlaneY, 'zdata', eyes.R.ListingPlaneZ);

        end

    end

    methods(Static) % Demo helpers
        function [worldPoints] = MakeWorldPoints(stimulusType, stimulusSizeCm, stimulusDistanceCm, stimulusTiltDeg, stimulusSlantDeg)
            persistent worldPointsNorm;
            persistent lastStimulusType;

            if ( isempty( lastStimulusType) || lastStimulusType ~= string(stimulusType) )
                lastStimulusType = stimulusType;
                sizeStimCm = 40;

                worldPoints = table();

                switch (stimulusType)
                    case 'RANDOMPLANE'
                        numDots = 200;
                        worldPoints.X = zeros(numDots, 1);
                        worldPoints.Y = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
                        worldPoints.Z = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
                    case 'GRIDPLANE'
                        [X,Y] = meshgrid(-(sizeStimCm/2):2:(sizeStimCm/2),-(sizeStimCm/2):2:(sizeStimCm/2));
                        worldPoints.X = zeros(size(X(:)));
                        worldPoints.Y = Y(:);
                        worldPoints.Z = X(:);

                    case 'HLINE'
                        [X,Y] = meshgrid(-(sizeStimCm/2):2:(sizeStimCm/2),(0)*ones(size(-(sizeStimCm/2):2:(sizeStimCm/2))));
                        worldPoints.X = zeros(size(X(:)));
                        worldPoints.Y = Y(:);
                        worldPoints.Z = X(:);
                    case 'VLINE'
                        [X,Y] = meshgrid((0)*ones(size(-(sizeStimCm/2):2:(sizeStimCm/2))),-(sizeStimCm/2):2:(sizeStimCm/2));
                        worldPoints.X = zeros(size(X(:)));
                        worldPoints.Y = Y(:);
                        worldPoints.Z = X(:);
                    case 'CROSS'
                        [X,Y] = meshgrid(-(sizeStimCm/2):2:(sizeStimCm/2),(0)*ones(size(-(sizeStimCm/2):2:(sizeStimCm/2))));
                        worldPoints.X = zeros(size(X(:)));
                        worldPoints.Y = X(:);
                        worldPoints.Z = Y(:);
                        [X,Y] = meshgrid((0)*ones(size(-(sizeStimCm/2):2:(sizeStimCm/2))),-(sizeStimCm/2):2:(sizeStimCm/2));
                        worldPoints2 = table();
                        worldPoints2.X = zeros(size(X(:)));
                        worldPoints2.Y = X(:);
                        worldPoints2.Z = Y(:);

                        worldPoints = cat(1,worldPoints, worldPoints2);
                end
                worldPointsNorm = worldPoints;
            else
                worldPoints = worldPointsNorm;
            end

            if (~isempty(worldPoints))
                % Resize stimulus
                worldPoints.X = zeros(size(worldPoints.Z));
                worldPoints.Y = worldPoints.Y*stimulusSizeCm/40;
                worldPoints.Z = worldPoints.Z*stimulusSizeCm/40;

                % Rotate by slant and tilt
                R = Geometry3D.AxisAngle2Mat([0 cosd(stimulusTiltDeg) sind(stimulusTiltDeg)],deg2rad(stimulusSlantDeg));
                worldPoints{:,:} = worldPoints{:,:}*R;

                % Displace by distance
                worldPoints.X = worldPoints.X + stimulusDistanceCm;
            end
        end
    end


    methods(Static) % DEMO MOTION FLOW

        function demoMotionFlow()

            app = InteractiveUI('Motion Flow Simulator',@(app) (Geometry3D.demoMotionFlowUpdate(app)), 0.2);
            % app.AddDropDown('Stimulus',      1,  ["Ground plane" "Point cloud" "Infinite"])
            app.AddSlider('Eye height',           1.5,[0    10])
            app.AddSlider('Fixation Azimuth',     0,  [-90 90])
            app.AddSlider('Fixation Elevation',   0,  [-90 90])
%             app.AddSlider('Ground plane slant',          0,  [-90  90])
            % app.AddSlider('Ground plane tilt',           0,  [0    90])
            app.AddSlider('Heading Speed',    0,  [0 10])
            app.AddSlider('Heading Azimuth',    0,  [-90 90])
            app.AddSlider('Heading Elevation',    0,  [-90 90])

            app.AddSlider('Head Angular velocity X',   0,  [-100 100])
            app.AddSlider('Head Angular velocity Y',   0,  [-100 100])
            app.AddSlider('Head Angular velocity Z',   0,  [-100 100])

            app.AddDropDown('Stabilize',          2,  ["YES" "NO"])
            app.AddSlider('Stabilization Gain',   1,  [-2 2])
            app.AddSlider('VORvsPursuitWeight',   1,  [0 1], ...
                'Relative contribution of stabilization and pursuit to stabilization, 0 - pure VOR, 1 - pure pursuit') % TODO: how to think about t-vor
            app.AddSlider('Listings law fraction',   0,  [0 1], ...
                'Tilt of angular velocity relative to Listing''s plane, 0 - angular velocity falls in listings plane, 1 - angular velocity is orthogonal to gaze direction')
            app.AddSlider('Eye Angular velocity X',   0,  [-100 100])
            app.AddSlider('Eye Angular velocity Y',   0,  [-100 100])
            app.AddSlider('Eye Angular velocity Z',   0,  [-100 100])
            app.AddDropDown('Coordinate system',      1,  [ "Fick" "Helmholtz" "Hess" "Harms"])
            % app.AddDropDown('View3D',             1,  ["Oblique" "TOP" "SIDE"])

            %
            % app.Data.Screen = struct();
            % app.Data.Screen.SizeCm = [30*16/9 30];
            % app.Data.Screen.ResPix = [1920 1080];
            % app.Data.Screen.Distance = 57;
            % app.Data.Screen.Slant = 0;
            %
            % app.Data.FixationSpot = struct();
            % app.Data.FixationSpot.X = 0;
            % app.Data.FixationSpot.Y = 0;
            % app.Data.FixationSpot.Z = 0;

            % Add average monocular motion flow and paralax flow

            app.Open();

        end

        function  h = demoMotionFlowInitPlots()

            % struct containing all the graphical objects
            h = struct();

            % create the figure
            screen_size = get(0,'ScreenSize');
            margin = floor(0.1*(screen_size(4)));
            h.figure = figure('color','w','position',floor([margin margin screen_size(3)*2.8/4 screen_size(4)*2/4 ]));

            plotHeadCentered        = [0.1    0.1    0.4    0.9];
            plotEyeCentered         = [0.6    0.1    0.4    0.9];
            plotEyeCenteredCoord    = [0.7    0.1    0.3    0.9];


            numberOfPoints = 1000;
            h.VisualDirections = Geometry3D.SampleVisualDirections(numberOfPoints,'Random');

            %-----------------------------
            % Head centered plot
            %-----------------------------
            h.plotH.ax = axes('OuterPosition', plotHeadCentered, 'nextplot','add');
            h.plotH.fixationSpot       = line(0,0, 'linestyle','none','marker','o','Color','r','LineWidth',2, 'markersize',10);
            h.plotH.headingDirection   = line(0,0, 'linestyle','none','marker','o','Color','g','LineWidth',2, 'markersize',6);
            h.plotH.motionFlow          = quiver(0,0,0,0,'autoscale','off', 'MaxHeadSize', 0.025);
            title({'Head centered'});

            axis equal; % Ensure equal scaling on both axes
            set(gca,"XLim",[-90 90],"YLim",[-90 90])
            grid off; axis off;

            % Define parameters for circles and radial lines
            R = 90;
            numCircles = 6;                   % Number of concentric circles
            numRadialLines = 12;              % Number of radial lines
            labelAngle = pi/8;                % Position labels at 45 degrees for clarity
            Geometry3D.DrawPolarGrid(R, numCircles, numRadialLines, labelAngle, 0 );




%             grid
%             ylim([-100 100]), zlim([-100 100]), xlim([-5 200]);
%             xlabel('X (cm)'), zlabel('Z (cm)'), ylabel('Y (cm)');
%             % view(-60,10)
%             title('3D world');
% 
%             % world points, shadow, and fixation
%             h.plot3D.wordlPoints        = line(0,0, 0,'linestyle','none','marker','o','Color','k');
%             h.plot3D.wordlPointsShadow  = line(0,0, 0,'linestyle','none','marker','o','Color',[0.8 0.8 0.8]);
%             h.plot3D.fixationSpot       = line(0,0, 0,'linestyle','none','marker','o','Color','r','LineWidth',2, 'markersize',20);
% 
%             % eyes
%             h.plot3D.LeftEye.Center     = plot3(0, 0, 0,'o','Color','b', 'markersize',10); % left eye fixation spot and right eye
%             h.plot3D.RightEye.Center    = plot3(0, 0, 0,'o','Color','r', 'markersize',10); % left eye fixation spot and right eye
%             h.plot3D.LeftEye.Lines      = line(0,0,0,'linestyle','-','Color','b'); 
%             h.plot3D.RightEye.Lines     = line(0,0,0,'linestyle','-','Color','r'); 
% 
%             % horopter surface
%             h.plot3D.horopter = surf(zeros(2),zeros(2),zeros(2),zeros(2));
%             originalCmap = colormap('turbo');
%             set(h.plot3D.horopter,'EdgeColor','none');
%             set(h.plot3D.horopter,'FaceAlpha',0.7);
%             set(h.plot3D.ax,'Projection','perspective')
%             set(h.plot3D.ax,'CameraPosition',[-100 -100 100],'cameratarget',[200 0 0])
% 
%             nColors = size(originalCmap, 1);
%             % Create a normalized index vector from 0 to 1
%             xi = linspace(0, 1, nColors);
%             xi2 = linspace(0, 1, nColors*4);
% 
%             % Apply a transformation to compact one side.
%             % Here, squaring the indices (xi.^2) compresses values near 0 more than near 1.
%             xi_transformed = xi2.^0.5;  % Adjust the exponent as needed
% 
%             % Interpolate to create the new colormap based on the transformed indices
%             newCmap = interp1(xi, originalCmap, xi_transformed);
%             colormap(newCmap);
%             set(gca,'clim',[0 50])
% 
%             % legend for 3d plot
%             h.plot3D.legend = legend([h.plot3D.LeftEye.Lines, h.plot3D.RightEye.Lines, h.plot3D.fixationSpot, h.plot3D.wordlPoints], {'Left eye','Right eye','Fixation point','Stimulus'},'Location','southwest','box','off');

            %-----------------------------
            % retina centered plot
            %-----------------------------
            
            h.plotR.ax = axes('OuterPosition', plotEyeCentered, 'nextplot','add');
            title({'Retina centered'});
        

            axis equal; % Ensure equal scaling on both axes
            set(gca,"XLim",[-90 90],"YLim",[-90 90])
            grid off; axis off;

            % Define parameters for circles and radial lines
            R = 90;
            numCircles = 6;                   % Number of concentric circles
            numRadialLines = 12;              % Number of radial lines
            labelAngle = pi/8;                % Position labels at 45 degrees for clarity
            Geometry3D.DrawPolarGrid(R, numCircles, numRadialLines, labelAngle, 0);


            h.plotR.fixationSpot       = line(0,0, 'linestyle','none','marker','o','Color','r','LineWidth',2, 'markersize',10);
            h.plotR.headingDirection   = line(0,0, 'linestyle','none','marker','o','Color','g','LineWidth',2, 'markersize',6);
            h.plotR.motionFlow          = quiver(0,0,0,0,'autoscale','off', 'MaxHeadSize', 0.025);

%             axis equal; % Ensure equal scaling on both axes
%             set(gca,"XLim",[-90 90],"YLim",[-90 90])
%             grid off; axis off;
%             set(gca,'clim',[0 50])
%             c = colorbar;
%             c.Label.String = 'Min Disparity (deg)';
% 
%             % Define parameters for circles and radial lines
%             R = 90;
%             numCircles = 6;                   % Number of concentric circles
%             numRadialLines = 12;              % Number of radial lines
%             labelAngle = pi/8;                % Position labels at 45 degrees for clarity
%             Geometry3D.DrawPolarGrid(R, numCircles, numRadialLines, labelAngle );
% 
%             h.plotRetina.horopter = surf(zeros(2),zeros(2),zeros(2),zeros(2));
%             colormap(newCmap);
%             h.plotRetina.horopter.EdgeColor = 'none';
%             h.plotRetina.horopter.FaceAlpha = 0.4;
% 
%             % stimulus in the left and right eye
%             h.plotRetina.leftEyePoints= plot(0, 0,'bo');
%             h.plotRetina.rightEyePoints = plot(0, 0,'ro');
% 
%             % fixation point
%             h.plotRetina.fixationSpot = plot(0,0,'ro','linewidth',3); % TODO: add left and right eye fixation point for missaligned eyes
%             legend([h.plotRetina.leftEyePoints, h.plotRetina.rightEyePoints , h.plotRetina.fixationSpot], {'Left eye','Right eye','Fixation point'},'Location','northeast','box','off')


            %-----------------------------
            % Plot retina coordinate system plot
            %-----------------------------

            % h.plotRC.ax = axes('OuterPosition', plotEyeCenteredCoord, 'nextplot','add');
%             set(gca, 'PlotBoxAspectRatio',[1 1 1])
%             grid
%             set(gca,'xlim',[-40 40],'ylim',[-40 40])
%             xlabel('Horizontal disparity (Helmtholz, deg, positive uncrossed)')
%             ylabel('Vertical disparity (Helmtholz, deg)')
% 
%             h.plotDisparity.disparities = quiver(0,0, 0, 0, 'AutoScale', "off");
% 
%             %-----------------------------
%             % haploscope plot
%             %-----------------------------
% 
%             h.plotHaploscope = [];
% 
%             h.plotHaploscope.axLeft = axes('OuterPosition', plotHaploscopeLeftPosition, 'nextplot','add');
%             h.plotHaploscope.leftPoints     = plot(0, 0, 'bo');
%             h.plotHaploscope.leftFixation   = plot(0, 0, 'ro', 'linewidth',2,'markersize',15);
%             grid
%             set(gca,'xlim',[0 1920],'ylim',[0 1080])
%             set(gca,'PlotBoxAspectRatio',[16 9 1])
%             xlabel({'Haploscope render left eye screen' '53.3cm x 30 cm, 1920x80, 57 cm dist'});
%             set(gca,'box','on')
% 
%             h.plotHaploscope.axRight = axes('OuterPosition',plotHaploscopeRightPosition, 'nextplot','add');
%             h.plotHaploscope.rightPoints  = plot(0, 0, 'ro');
%             h.plotHaploscope.rightFixation  = plot(0, 0, 'ro', 'linewidth',2,'markersize',15);
%             grid
%             set(gca,'xlim',[0 1920],'ylim',[0 1080])
%             set(gca,'yticklabels',[])
%             set(gca,'PlotBoxAspectRatio',[16 9 1])
%             xlabel({'Haploscope render right eye screen' '53.3cm x 30 cm, 1920x80, 57 cm dist'});
%             set(gca,'box','on')
% 
%             %-----------------------------
%             % Listing planes plot
%             %-----------------------------
% 
%             % eyes
%             h.plotListings.ax = axes('OuterPosition',plotListingsPosition, 'nextplot','add');
%             view(h.plotListings.ax, -40 ,  40)
%             axis equal
%             grid
%             set(h.plotListings.ax,'xlim',[-15 15],'ylim',[-15 15],'zlim',[-15 15])
%             title('Listing''s planes')
% 
%             % eyes
%             h.plotListings.LeftEye.Center     = plot3(0, 0, 0,'o','Color','b', 'markersize',10); % left eye fixation spot and right eye
%             h.plotListings.RightEye.Center    = plot3(0, 0, 0,'o','Color','r', 'markersize',10); % left eye fixation spot and right eye
%             h.plotListings.LeftEye.Lines      = line(0,0,0,'linestyle','-','Color','b'); 
%             h.plotListings.RightEye.Lines     = line(0,0,0,'linestyle','-','Color','r'); 
% 
%             % Listing planes surface
%             h.plotListings.LeftPlane    = surf(zeros(2),zeros(2),zeros(2),zeros(2),'facecolor','white','edgecolor','b');
%             h.plotListings.RightPlane   = surf(zeros(2),zeros(2),zeros(2),zeros(2),'facecolor','white','edgecolor','r');

        end

        function demoMotionFlowUpdate(app)
                
%                Stimulus: "Ground plane"
%               EyeHeight: 150
%         FixationAzimuth: 0
%       FixationElevation: 0
%            HeadingSpeed: 0
%          HeadingAzimuth: 0
%        HeadingElevation: 0
%               Stabilize: "NO"
%       StabilizationGain: 1
%     ListingsLawFraction: 0
%        AngularVelocityX: 0
%        AngularVelocityY: 0
%        AngularVelocityZ: 0
%        CoordinateSystem: "Fick"


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update graphics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if ( ~isfield(app.Data, "h") || ~isfield(app.Data.h, "figure") || ~isvalid(app.Data.h.figure))
                % If figure does not exist create it with all the plots and
                % handles to them
                app.Data.h = Geometry3D.demoMotionFlowInitPlots();
            end

            coordinateSystem = app.Values.CoordinateSystem;

            % heading velocity in head reference
            [x,y,z] = Geometry3D.AzElToSphere( deg2rad(app.Values.HeadingAzimuth), deg2rad(app.Values.HeadingElevation), coordinateSystem);
            headingDirection = [x,y,z]';
            headingVelocity = app.Values.HeadingSpeed*headingDirection;

            [x,y,z] = Geometry3D.AzElToSphere( deg2rad(app.Values.FixationAzimuth), deg2rad(app.Values.FixationElevation), coordinateSystem);
            gazeDirection = [x,y,z]';
            eyePositionRotMat = Geometry3D.LookAtListingsSimple(gazeDirection);
            eyeLocation = [0, 0, app.Values.EyeHeight];
            eyeLinearVelocity = eyePositionRotMat'*headingVelocity;

            headAngularVelocity = deg2rad([app.Values.HeadAngularVelocityX app.Values.HeadAngularVelocityY app.Values.HeadAngularVelocityZ]');
            eyeheadAngularVelocity = eyePositionRotMat'*headAngularVelocity;
    
            foveaDirection = [1 0 0];
            fovealVelocity = Geometry3D.CalculateMotionField(foveaDirection, eyeheadAngularVelocity, eyeLinearVelocity, eyeLocation, eyePositionRotMat);
            [JacAngular] = Geometry3D.CalculateMotionJacobianFields(foveaDirection);

            gain = app.Values.StabilizationGain;
            if app.Values.Stabilize == "YES"
                eyeAngularVelocity = -gain*JacAngular'*fovealVelocity'; 

                normal = gazeDirection * (1-app.Values.ListingsLawFraction) + [1 0 0]' * (app.Values.ListingsLawFraction);
                wx = -(eyeAngularVelocity(2)*normal(2) + eyeAngularVelocity(3)*normal(3))/normal(1);
                eyeAngularVelocity(1) = wx;
                eyeAngularVelocity = eyeAngularVelocity*(app.Values.VORvsPursuitWeight) - eyeheadAngularVelocity*(1-app.Values.VORvsPursuitWeight);
            else
                eyeAngularVelocity = [0 0 0]';
            end




            eyeAngularVelocity = eyeAngularVelocity + deg2rad([app.Values.EyeAngularVelocityX app.Values.EyeAngularVelocityY app.Values.EyeAngularVelocityZ]');

            
            visualDirections = app.Data.h.VisualDirections;


            % 1 - head reference motion field (no eye rotation)
            % 2 - eye reference motion field without eye velocity
            % 3 - eye reference motion field with eye velocity
            motionFieldLinearHeadReference = Geometry3D.CalculateMotionField(visualDirections, headAngularVelocity, headingVelocity, eyeLocation);
            motionFieldLinearEyeReference = Geometry3D.CalculateMotionField(visualDirections, eyeheadAngularVelocity, eyeLinearVelocity, eyeLocation, eyePositionRotMat);
            motionFieldTotalEyeRefefence = Geometry3D.CalculateMotionField(visualDirections, eyeheadAngularVelocity+eyeAngularVelocity, eyeLinearVelocity, eyeLocation, eyePositionRotMat);


            
            gazeDirectionProjected = rad2deg(Geometry3D.AzimuthalEquidistantProjectionZY(gazeDirection'));
            headingDirectionProjected = rad2deg(Geometry3D.AzimuthalEquidistantProjectionZY(headingDirection'));
            visualDirectionsProjected = rad2deg(Geometry3D.AzimuthalEquidistantProjectionZY(visualDirections));


            h = app.Data.h;

            set(h.plotH.fixationSpot, 'xdata', gazeDirectionProjected(2), 'ydata', gazeDirectionProjected(3) );
            set(h.plotH.headingDirection, 'xdat', headingDirectionProjected(2), 'ydata', headingDirectionProjected(3) );
            set(h.plotH.motionFlow, 'xdata', visualDirectionsProjected(:,2), 'ydata', visualDirectionsProjected(:,3), 'udata', rad2deg( motionFieldLinearHeadReference(:,1)), 'vdata', rad2deg(motionFieldLinearHeadReference(:,2)));


            gazeDirectionProjected = rad2deg(Geometry3D.AzimuthalEquidistantProjectionZY([1 0 0]));
            headingDirectionProjected = rad2deg(Geometry3D.AzimuthalEquidistantProjectionZY(headingDirection'*eyePositionRotMat));
            visualDirectionsProjected = rad2deg(Geometry3D.AzimuthalEquidistantProjectionZY(visualDirections));

            set(h.plotR.fixationSpot, 'xdata', gazeDirectionProjected(2), 'ydata', gazeDirectionProjected(3) );
            set(h.plotR.headingDirection, 'xdat', headingDirectionProjected(2), 'ydata', headingDirectionProjected(3) );
            set(h.plotR.motionFlow, 'xdata', visualDirectionsProjected(:,2), 'ydata', visualDirectionsProjected(:,3), 'udata', rad2deg(motionFieldTotalEyeRefefence(:,1)), 'vdata', rad2deg(motionFieldTotalEyeRefefence(:,2)));

% 
% 
%             Values = app.Values;
% 
%             if (~isfield(app.Data,"stimulus"))
%                 app.Data.stimulus = "NONE";
%             end
%             if (~isfield(app.Data,"View3D"))
%                 app.Data.View3D = "NONE";
%             end
% 
%             %%
%             sizeStimCm = 40;
% 
%             %             ["Ground plane" "Cloud"])
% 
%             if (app.Data.stimulus ~= string(Values.Stimulus) )
%                 % Make some dots in world coordinates (XYZ)
%                 %
%                 % Make a table to start putting everything together
%                 switch (Values.Stimulus)
%                     case 'Cloud'
%                         numDots = 200;
%                         worldPoints = table();
%                         worldPoints.X = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
%                         worldPoints.Y = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
%                         worldPoints.Z = zeros(numDots, 1);
%                     case 'Ground plane'
%                         % TODO: need to update the update funciton so the Z is not overwritten by a
%                         % constant distance. Just shift things.
%                         [X,Y] = meshgrid([-(sizeStimCm/2):2:(sizeStimCm/2)],[-(sizeStimCm/2):2:(sizeStimCm/2)]);
%                         worldPoints = table();
%                         worldPoints.X = X(:);
%                         worldPoints.Y = Y(:);
%                         worldPoints.Z = zeros(size(worldPoints.X));
%                 end
%                 app.Data.worldPoints = worldPoints;
%                 app.Data.stimulus = Values.Stimulus;
%             else
%                 worldPoints = app.Data.worldPoints;
%             end
% 
% 
% 
%             %% Update the points, eyes and screen acordint to the sliders
%             app.Data.FixationSpot.X = Values.FixationX;
%             app.Data.FixationSpot.Y = Values.FixationY;
%             app.Data.FixationSpot.Z = Values.FixationDistance;
%             app.Data.Screen.Slant = Values.ScreenSlant;
% 
%             leftEyeScreen = Geometry3D.MakeScreen(app.Data.Screen.SizeCm, app.Data.Screen.ResPix, app.Data.Screen.Distance, app.Data.Screen.Slant);
%             rightEyeScreen = Geometry3D.MakeScreen(app.Data.Screen.SizeCm, app.Data.Screen.ResPix, app.Data.Screen.Distance, app.Data.Screen.Slant);
% 
%             eyes = Geometry3D.MakeEyes(Values.IPDMm/10, app.Data.FixationSpot, Values.TorsionVersion, Values.TorsionVergence);
% 
%             % Rotate the world points to apply the slant (rotates around
%             % 0,0,0)
%             worldPoints.X = worldPoints.X*Values.StimulusScale;
%             worldPoints.Y = worldPoints.Y*Values.StimulusScale;
%             worldPoints.Z = zeros(size(worldPoints.Z));
% 
%             % Rotate by slant and tilt
%             R = Geometry3D.Quat2RotMat(Geometry3D.AxisAngle2Quat([cosd(Values.PlaneTilt) sind(Values.PlaneTilt) 0],deg2rad(Values.PlaneSlant)));
%             worldPoints{:,:} = (R*worldPoints{:,:}')';
% 
%             % Displace by distance
%             worldPoints.Z = worldPoints.Z + Values.StimulusDistance;
% 
%             % Add the fixation spot to the world points to conver to eye
%             % and screen points.
%             worldPoints{end+1,:} = [app.Data.FixationSpot.X app.Data.FixationSpot.Y app.Data.FixationSpot.Z];
% 
%             %% Get the points in eye and screen coordinates
%             eyePoints = Geometry3D.Points3DToEyes(worldPoints, eyes);
%             screenPoints = Geometry3D.PointsEyesToScreen(eyes, eyePoints, leftEyeScreen, rightEyeScreen);
% 
% 
%             % update 3D plot
%             lxfar = eyes.L.X - 10000*eyes.L.RotMat(2,1);
%             lyfar = eyes.L.Center(2) - 10000*eyes.L.RotMat(3,1);
%             lzfar = eyes.L.Center(3) + 10000*eyes.L.RotMat(1,1);
% 
%             rxfar = eyes.R.Center(1) - 10000*eyes.R.RotMat(2,1);
%             ryfar = eyes.R.Center(2) - 10000*eyes.R.RotMat(3,1);
%             rzfar = eyes.R.Center(3) + 10000*eyes.R.RotMat(1,1);
% 
%             lrx = eyes.L.X - 10*eyes.L.RotMat(2,2);
%             llx = eyes.L.X + 10*eyes.L.RotMat(2,2);
%             lry = eyes.L.Center(2) - 10*eyes.L.RotMat(3,2);
%             lly = eyes.L.Center(2) + 10*eyes.L.RotMat(3,2);
%             lrz = eyes.L.Center(3) + 10*eyes.L.RotMat(1,2);
%             llz = eyes.L.Center(3) - 10*eyes.L.RotMat(1,2);
% 
%             rrx = eyes.R.Center(1) - 10*eyes.R.RotMat(2,2);
%             rlx = eyes.R.Center(1) + 10*eyes.R.RotMat(2,2);
%             rry = eyes.R.Center(2) - 10*eyes.R.RotMat(3,2);
%             rly = eyes.R.Center(2) + 10*eyes.R.RotMat(3,2);
%             rrz = eyes.R.Center(3) + 10*eyes.R.RotMat(1,2);
%             rlz = eyes.R.Center(3) - 10*eyes.R.RotMat(1,2);
% 
%             rbx = eyes.R.Center(1) - 10*eyes.R.RotMat(2,3);
%             rtx = eyes.R.Center(1) + 10*eyes.R.RotMat(2,3);
%             rby = eyes.R.Center(2) - 10*eyes.R.RotMat(3,3);
%             rty = eyes.R.Center(2) + 10*eyes.R.RotMat(3,3);
%             rbz = eyes.R.Center(3) + 10*eyes.R.RotMat(1,3);
%             rtz = eyes.R.Center(3) - 10*eyes.R.RotMat(1,3);
% 
%             lbx = eyes.L.X - 10*eyes.L.RotMat(2,3);
%             ltx = eyes.L.X + 10*eyes.L.RotMat(2,3);
%             lby = eyes.L.Center(2) - 10*eyes.L.RotMat(3,3);
%             lty = eyes.L.Center(2) + 10*eyes.L.RotMat(3,3);
%             lbz = eyes.L.Center(3) + 10*eyes.L.RotMat(1,3);
%             ltz = eyes.L.Center(3) - 10*eyes.L.RotMat(1,3);
% 
%             set(app.Data.hfix, 'xdata', Values.FixationX);
%             set(app.Data.hfix, 'ydata', Values.FixationDistance);
%             set(app.Data.hfix, 'zdata', Values.FixationY);
%             set(app.Data.heyes.l, ...
%                 'xdata', [eyes.L.X lxfar ...
%                 nan lbx ltx nan lrx llx], ...
%                 'ydata', [eyes.L.Center(3) lzfar  ...
%                 nan lbz ltz nan lrz llz], ...
%                 'zdata', [eyes.L.Center(2) lyfar ...
%                 nan lby lty nan lry lly]);
%             set(app.Data.heyes.r, ...
%                 'xdata', [ eyes.R.Center(1) rxfar ...
%                 nan rbx rtx nan rrx rlx ], ...
%                 'ydata', [ eyes.R.Center(3) rzfar ...
%                 nan rbz rtz nan rrz rlz ], ...
%                 'zdata', [ eyes.R.Center(2) ryfar ...
%                 nan rby rty nan rry rly]);
%             set(app.Data.heyes.c(1), 'xdata', [eyes.L.X]);
%             set(app.Data.heyes.c(2), 'xdata', [eyes.R.Center(1)]);
% 
%             set(app.Data.hpoints, 'xdata', worldPoints.X, 'ydata',worldPoints.Z, 'zdata', worldPoints.Y);
%             set(app.Data.hspoints, 'xdata', worldPoints.X, 'ydata',worldPoints.Z, 'zdata', app.Data.hspoints.Parent.XLim(1)*ones(size(worldPoints.Y)));
% 
%             ax3d = app.Data.hfix.Parent;
%             switch(app.Values.View3D)
%                 case 'Oblique'
%                     view(ax3d, -60,10)
%                 case 'TOP'
%                     view(ax3d, -90,90)
%                 case 'SIDE'
%                     view(ax3d, -90,0)
%             end
% 
%             % update retina plot
%             ra = atan2( tand(eyePoints.RV) , tand(eyePoints.RH)./(abs(cosd(eyePoints.RV))));
%             re = atand(sqrt(tand(eyePoints.RV).^2 + (tand(eyePoints.RH)./(abs(cosd(eyePoints.RV)))).^2));
% 
%             set(app.Data.hLRpoints.R, 'ThetaData', ra, 'RData', re);
%             ra = atan2( tand(eyePoints.LV) , tand(eyePoints.LH)./abs(cosd(eyePoints.LV)));
%             re = atand(sqrt(tand(eyePoints.LV).^2 + (tand(eyePoints.LH)./(abs(cosd(eyePoints.LV)))).^2));
%             set(app.Data.hLRpoints.L, 'ThetaData', ra, 'RData', re);
% 
%             % update disparity plot
%             set(app.Data.hdisparity, ...
%                 'xdata',eyePoints.RH, 'ydata', eyePoints.RV, ...
%                 'udata', eyePoints.LH-eyePoints.RH, 'vdata', eyePoints.LV-eyePoints.RV);
% 
%             % update screen plots
%             set(app.Data.hscreen.LPoints, 'xdata', screenPoints.LX(1:end-1), 'ydata',screenPoints.LY(1:end-1));
%             set(app.Data.hscreen.RPoints, 'xdata', screenPoints.RX(1:end-1), 'ydata',screenPoints.RY(1:end-1));
%             set(app.Data.hscreen.LFP, 'xdata', screenPoints.LX(end), 'ydata', screenPoints.LY(end));
%             set(app.Data.hscreen.RFP, 'xdata', screenPoints.RX(end), 'ydata', screenPoints.RY(end));

        end

    end

    methods(Static) % DEMO LISTINGS LAW

        function demoListingsLaw()
            app = InteractiveUI('Listing''s law demo',@(app) (Geometry3D.demoListingsLawUpdate(app)), .2);

            app.AddDropDown('Coordinate system',      1,  [ "Fick" "Helmholtz" "Hess" "Harms" "ImagePlane"])
            app.AddSlider('Eye Position Azimuth',     0, [-60 60], ...
                'Eye position horizontal motion in whichever coordinate system is selected.')
            app.AddSlider('Eye position Elevation',    0, [-60 60], ...
                'Eye position vertical motion in whichever coordinate system is selected.')
            app.AddSlider('Torsion',0, [-60 60], ...
                'Eye torsion. If following Listing''s law it just gets added on top.')

            app.AddSlider('Listing''s plane pitch',   0,  [-90 90])
            app.AddSlider('Listing''s plane yaw',     0,  [-90 90])

            app.AddDropDown('Follow Listing''s Law?',           2,  ["YES" "NO"])

            app.AddSlider('World fixed image patch size',   0,  [0 10])
            app.AddSlider('World fixed patch azimuth',        0,  [-60 60])
            app.AddSlider('World fixed patch elevation',      5,  [-60 60])
            app.AddSlider('Retina fixed image patch size',   0,  [0 10])
            app.AddSlider('Retina fixed patch azimuth',       5,  [-60 60])
            app.AddSlider('Retina fixed patch elevation',     0,  [-60 60])

            app.Open();
        end

        function [f, h] = demoListingsLawInitPlots(app)
            scr_siz = get(0,'ScreenSize');
            margin = floor(0.1*(scr_siz(4)));
            f = figure('color','w','position',floor([...
                margin...
                margin...
                scr_siz(3)*2/4 ...
                scr_siz(4)*2/4 ...
                ]));

            % set(gcf, 'Renderer', 'OpenGL');
             % shading interp, material shiny, lighting phong, lightangle(0, 55);

            % Define the radius of the sphere
            R = 1;
            % Define the resolution of the sphere
            res = 50;
            % Define the transparency of the sphere
            alpha = 1;

            %
            % Create the sphere
            [Z,X,Y] = sphere(res);

            % Plot the sphere
            h.sphere = mesh(R*X,R*Y,R*Z,'FaceAlpha', alpha,'facecolor',[0.8 0.8 0.8],'edgecolor','none');
            hold
            h.sphere2 = mesh(R*X,R*Y,R*Z,'FaceAlpha', 0,'facecolor',[0.8 0.8 0.8],'edgecolor',0.2*[0.8 0.8 0.8]);

            h.screen = mesh(R*X,R*Y,R*Z,'FaceAlpha', 0,'facecolor',[0.8 0.8 0.8],'edgecolor',0.2*[0.8 0.8 0.8]);
            h.listingsPlane =  mesh(R*X,R*Y,R*Z,'FaceAlpha', 0,'facecolor',[0.8 0.8 0.8],'edgecolor',[0.8 0.8 0.8]);
            h.rotationVector = line(0,0,0,'linewidth',4,'color',[0.3 0.6 0.3], 'linestyle','-.');
            h.rotationVectorFromPrimary = line(0,0,0,'linewidth',4,'color',[0.6 0.3 0.3], 'linestyle','--');
            axis equal;

    
            set(gca,'OuterPosition', [-0.72   0   1.7     1])
            h.ax = gca;

            h.frame(1) = line(0,0,0,'linewidth',2,'color','k');
            h.frame(2) = line(0,0,0,'linewidth',2,'color','g');
            % h.frame(3) = line(0,0,0,'linewidth',2,'color','k');


            h.Vmeridian = line(0,0,0,'linewidth',2,'color','r');
            h.Hmeridian = line(0,0,0,'linewidth',2,'color','b');

            h.rVmeridianScreen = line(0,0,0,'linewidth',2,'color','r');
            h.rHmeridianScreen = line(0,0,0,'linewidth',2,'color','b');          
            h.wVmeridianScreen = line(0,0,0,'linewidth',2,'color','g');
            h.wHmeridianScreen = line(0,0,0,'linewidth',2,'color','c');            

            set(gca,'Projection','perspective')
            set(gca,'CameraPosition',[-1 10 1],'cameratarget',[1 0 0])
            set(gca,'xlim',[-2 3],'ylim',[-4 4],'zlim',[-4 4])
            set(gca,'XTick',[],'YTick',[],'ZTick',[])
            %set(gca,'visible','off')
            grid off
            set(gca,'YDir','reverse')

            h.image3Dplot.WorldFixedPath = line(0,0,0,'linewidth',2,'color','r');
            h.image3Dplot.RetinaFixedPath = line(0,0,0,'markersize',6,'linewidth',2,'color','b');


            xlabel(h.ax,'x')
            ylabel(h.ax,'y')
            zlabel(h.ax,'z')

            legend([ h.listingsPlane,  h.rotationVector, h.rotationVectorFromPrimary,],...
                {'Listing''s Plane', 'Rotation axis from straight ahead', 'Rotation axis from Primary (listing''s) position'}, ...
                'Location','none','box','off','fontsize',12, 'Position',[0.12 0.9 0.08 0.1]);


            % screen plot

            screenD = 3;

            h.imageplot.ax = axes('OuterPosition', [0.3 0 0.35 1], 'nextplot','add');
            axis equal
            h.imageplot.screen = Geometry3D.Plot2DMeshGrid(Y./X*screenD,Z./X*screenD);
            set(h.imageplot.screen,'color',0.2*[0.8 0.8 0.8])
            h.imageplot.rVmeridianScreen = line(0,0,'linewidth',2,'color','r');
            h.imageplot.rHmeridianScreen = line(0,0,'linewidth',2,'color','b');
            h.imageplot.wVmeridianScreen = line(0,0,'linewidth',2,'color','g');
            h.imageplot.wHmeridianScreen = line(0,0,'linewidth',2,'color','c');
            h.imageplot.PrimaryPosition = line(0,0,'linestyle','none','marker','o','markersize',10,'linewidth',2,'color','g');
            h.imageplot.FixationPosition = line(0,0,'linestyle','none','marker','o','markersize',6,'linewidth',2,'color','k');
            h.imageplot.WorldFixedPath = line(0,0,'linewidth',2,'color','r');
            h.imageplot.RetinaFixedPath = line(0,0,'linewidth',2,'color','b');
            h.imageplot.WorldFixedPathApprox = line(0,0,'linewidth',1,'color','r');
            h.imageplot.RetinaFixedPathApprox = line(0,0,'linewidth',1,'color','b');
            set(gca,'xlim',deg2rad([-50 50])*screenD,'ylim',deg2rad([-50 50])*screenD)
            set(gca,'XTick',[],'YTick',[])
            title('Projection onto flat screen');

            legend([ h.imageplot.FixationPosition,  h.imageplot.PrimaryPosition, h.imageplot.wHmeridianScreen, h.imageplot.wVmeridianScreen , h.imageplot.rHmeridianScreen, h.imageplot.rVmeridianScreen ],...
                {'Fixation position', 'Primary (listing''s) position', 'Screen Horizontal meridian', 'Screen Vertical meridian', 'Retinal Horizontal meridian', 'Retinal Vertical meridian' }, ...
                'Location','none','box','off','fontsize',12, 'Position',[0.4 0.82 0.14 0.16]);

            h.imageplot.SVVtext = text(-2.2, -3.4,'hi');

            % retina plot

            h.imageplot2.ax = axes('OuterPosition', [0.65 0 0.35 1], 'nextplot','add');
            h.imageplot2.rVmeridianScreen = line(0,0,'linewidth',2,'color','r');
            h.imageplot2.rHmeridianScreen = line(0,0,'linewidth',2,'color','b');
            h.imageplot2.wVmeridianScreen = line(0,0,'linewidth',2,'color','g');
            h.imageplot2.wHmeridianScreen = line(0,0,'linewidth',2,'color','c');

            h.imageplot2.WorldFixedPath = line(0,0,'linewidth',2,'color','r');
            h.imageplot2.RetinaFixedPath = line(0,0,'markersize',6,'linewidth',2,'color','b');
            h.imageplot2.WorldFixedPathApprox = line(0,0,'linewidth',1,'color','r');
            h.imageplot2.RetinaFixedPathApprox = line(0,0,'linewidth',1,'color','b');
            h.imageplot2.FixationPosition = line(0,0,'linestyle','none','marker','o','markersize',6,'linewidth',2,'color','k');

            h.imageplot2.ErrorText = text(-40, -60,'hi');

            h.WorldFixedPathData = [[(-10:10)';10*ones(21,1);(10:-1:-10)';-10*ones(21,1)] [10*ones(21,1);(10:-1:-10)';-10*ones(21,1);(-10:10)' ]  ]/10;
            h.RetinaFixedPathData = [[(-10:10)';10*ones(21,1);(10:-1:-10)';-10*ones(21,1)] [10*ones(21,1);(10:-1:-10)';-10*ones(21,1);(-10:10)' ]  ]/10;


            legend([ h.imageplot2.RetinaFixedPath,  h.imageplot2.RetinaFixedPathApprox, ],...
                {'Retina fixed patch' 'Retina fixed patch approximation'}, ...
                'Location','none','box','off','fontsize',12, 'Position',[0.8 0.82 0.14 0.12]);
            
            % h.imageplot2.PrimaryPosition = line(0,0,'linestyle','none','marker','o','markersize',10,'linewidth',2,'color','g');
            % set(gca,'XTick',[],'YTick',[])
            axis equal
            set(gca,'xlim',4*[-10 10],'ylim',4*[-10 10])
            grid
            title('Projection onto retina coordinate system');
            xlabel('Azimuth (deg)'), ylabel('Elevation (deg)')

            % legend([  h.imageplot2.points2,  h.imageplot2.points3], {'Points shifted by full rotation', 'Points shifted by only translation (small angle aprox.)' },'Location','southoutside','box','off','fontsize',12);

        end

        function demoListingsLawUpdate(app)

            if ( ~isfield(app.Data, "f") || ~isvalid(app.Data.f))
                % If figure does not exist create it with all the plots and
                % handles to them
                [f, h] = Geometry3D.demoListingsLawInitPlots(app);
                app.Data.f = f;
                app.Data.h = h;
            end


            % screen distance
            screenD = 3;

            Values = app.Values;
            hvt = [Values.EyePositionAzimuth Values.EyePositionElevation Values.Torsion];


            % points for the grid
            steps = -85:5:85;
            rpoints = -85:5:85;
            [azgrid, elgrid] = meshgrid(steps,rpoints);

            [xgrid, ygrid, zgrid] = Geometry3D.AzElToSphere(deg2rad(azgrid), deg2rad(elgrid), app.Values.CoordinateSystem);
            [fx, fy, fz] = Geometry3D.AzElToSphere(deg2rad(hvt(1)), deg2rad(hvt(2)), app.Values.CoordinateSystem);
            [px,py,pz] = Geometry3D.AzElToSphere(deg2rad(Values.Listing_sPlaneYaw), deg2rad(Values.Listing_sPlanePitch), app.Values.CoordinateSystem);

            n = [px,py,pz]';
            v = [fx fy fz]';

            % if we are following listings law we just use the azimuth and
            % elevation to decide the 3D rotationa and then add the torsion
            % as rotation off of listings plane. 
            if Values.FollowListing_sLaw_ == "YES"
                [R,~, RprimToV] = Geometry3D.LookAtListings(v,n);
                RT = Geometry3D.RotX(-deg2rad(hvt(3)));
                R = R*RT;
                RprimToV = RprimToV*RT;
            else            
                R = Geometry3D.EulerAngles2RotMat(deg2rad(hvt).*[1 -1 -1], app.Values.CoordinateSystem);
                RprimToV = R; % TODO:
            end

            % calculate the tilt of SVV and SVH
            % R is your 3×3 rotation matrix
            v = R * [1;0;0];      % rotated LOS
            u = R * [0;1;0];      % rotated vertical meridian
            n2 = cross(v,u);       % normal to P
            L = cross(n2, [1;0;0]);% intersection with yz‐plane
            SVH = -atan2d(L(3), L(2));   % in radians, ±

            % R is your 3×3 rotation matrix
            v = R * [1;0;0];      % rotated LOS
            u = R * [0;0;1];      % rotated vertical meridian
            n2 = cross(v,u);       % normal to P
            L = cross(n2, [1;0;0]);% intersection with yz‐plane
            SVV = atan2d(L(2), L(3));   % in radians, ±

            set(app.Data.h.imageplot.SVVtext,'string',sprintf('Tilt of vertical meridian = %0.1f deg \nTilt of horizontal meridian = %0.1f deg \n',SVV,SVH))

            % world and retinal meridians in world and retinal coordinagtes

            % points for meridians in the retina
            theta = (-85:1:85)';
            rinrVmeridianPoints = [cosd(theta) zeros(size(theta)) sind(theta) ];
            rinrHmeridianPoints = [cosd(theta) sind(theta) zeros(size(theta)) ];
            % points for merdians in world
            theta = (-85:1:85)'*4;
            winwVmeridianPoints = screenD*[ones(size(theta)) zeros(size(theta)) deg2rad(theta) ];
            winwHmeridianPoints = screenD*[ones(size(theta)) deg2rad(theta)  zeros(size(theta))];

            % rotate the points in the world to retinal coordinates
            winrVmeridianPoints = normalize(winwVmeridianPoints*R', 2, 'norm');
            winrHmeridianPoints = normalize(winwHmeridianPoints*R', 2, 'norm');
            % rotate the points in the retina to world coordiantes and
            % project to screen
            rinwVmeridianPoints = rinrVmeridianPoints*R';
            rinwHmeridianPoints = rinrHmeridianPoints*R';
            rinsVmeridianPoints = rinwVmeridianPoints./rinwVmeridianPoints(:,1)*screenD;
            rinsHmeridianPoints = rinwHmeridianPoints./rinwHmeridianPoints(:,1)*screenD;





            wpatch3D = app.Values.WorldFixedImagePatchSize*app.Data.h.WorldFixedPathData;
            wpatch3D(:,1) = wpatch3D(:,1) + app.Values.WorldFixedPatchAzimuth;
            wpatch3D(:,2) = wpatch3D(:,2) + app.Values.WorldFixedPatchElevation;
            wpatch3D = deg2rad(wpatch3D)*screenD;
            wpatch3D(:,3) = screenD;
            wpatch = wpatch3D(:,[1 2]);
            wpatch3D = wpatch3D(:,[3 1 2]);

            winrpatch3D = (R'*wpatch3D')';



            rpatch3D = app.Values.RetinaFixedImagePatchSize*app.Data.h.WorldFixedPathData;
            rpatch3D(:,1) = rpatch3D(:,1) + app.Values.RetinaFixedPatchAzimuth;
            rpatch3D(:,2) = rpatch3D(:,2) + app.Values.RetinaFixedPatchElevation;
            rpatch3D = deg2rad(rpatch3D)*screenD;
            rpatch3D(:,3) = screenD;
            rpatch3D = rpatch3D(:,[3 1 2]);
            rinrrpatch3D = rpatch3D;
            rpatch3D = (R*rpatch3D')';


            rpatch = rpatch3D(:,[2 3])./rpatch3D(:,1)*screenD;
            rpatch3D = rpatch3D(:,[1 2 3])./rpatch3D(:,1)*screenD;

            rinwpatchAprox = app.Values.RetinaFixedImagePatchSize*app.Data.h.WorldFixedPathData;
            calibrationAngle = 20;
            calibrationFactor = calibrationAngle/tand(calibrationAngle);
            % rinwpatchAprox(:,1) = rinwpatchAprox(:,1) + calibrationFactor*tand(app.Values.RetinaFixedPatchAzimuth) + calibrationFactor*tand(hvt(1));
            % rinwpatchAprox(:,2) = rinwpatchAprox(:,2) + calibrationFactor*tand(app.Values.RetinaFixedPatchElevation) + calibrationFactor*tand(hvt(2));
            rinwpatchAprox(:,1) = rinwpatchAprox(:,1) + calibrationFactor*tand(app.Values.RetinaFixedPatchAzimuth) + rad2deg(fy/fx);
            rinwpatchAprox(:,2) = rinwpatchAprox(:,2) + calibrationFactor*tand(app.Values.RetinaFixedPatchElevation) + rad2deg(fz/fx);
            rinwinrpatchAprox = deg2rad(rinwpatchAprox);

            rinwpatchAprox = deg2rad(rinwpatchAprox)*screenD;
            rinwinrpatchAprox(:,3) = 1;
            rinwinrpatchAprox(:,1) = rinwinrpatchAprox(:,1);
            rinwinrpatchAprox(:,2) = rinwinrpatchAprox(:,2);
            rinwinrpatchAprox = rinwinrpatchAprox(:,[3 1 2]);
            rinwinrpatchAprox = (R'*rinwinrpatchAprox')';


            winrpatchAprox = app.Values.WorldFixedImagePatchSize*app.Data.h.WorldFixedPathData;
            winrpatchAprox(:,1) = winrpatchAprox(:,1) + app.Values.WorldFixedPatchAzimuth - hvt(1);
            winrpatchAprox(:,2) = winrpatchAprox(:,2) + app.Values.WorldFixedPatchElevation - hvt(2);
            winrpatchAprox = deg2rad(winrpatchAprox);


            % Get azimuth and elevation for retina plot
            [vmaz, vmel] = Geometry3D.SphereToAzEl(winrVmeridianPoints(:,1), -winrVmeridianPoints(:,2), -winrVmeridianPoints(:,3), app.Values.CoordinateSystem);
            [hmaz, hmel] = Geometry3D.SphereToAzEl(winrHmeridianPoints(:,1), -winrHmeridianPoints(:,2), -winrHmeridianPoints(:,3), app.Values.CoordinateSystem);
            [wpatchaz, wpatchel] = Geometry3D.SphereToAzEl(winrpatch3D(:,1), -winrpatch3D(:,2), -winrpatch3D(:,3), app.Values.CoordinateSystem);
            [rpatchaz, rpatchel] = Geometry3D.SphereToAzEl(rinrrpatch3D(:,1), -rinrrpatch3D(:,2), -rinrrpatch3D(:,3), app.Values.CoordinateSystem);
            [rpatchazApprox, rpatchelApprox] = Geometry3D.SphereToAzEl(rinwinrpatchAprox(:,1), -rinwinrpatchAprox(:,2), -rinwinrpatchAprox(:,3), app.Values.CoordinateSystem);




            % update 3D plot

            set(app.Data.h.sphere2, 'xdata', xgrid, 'ydata',ygrid, 'zdata',zgrid)
            set(app.Data.h.screen, 'xdata', xgrid./xgrid*screenD, 'ydata',ygrid./xgrid*screenD, 'zdata',zgrid./xgrid*screenD)

            set(app.Data.h.frame(1), 'xdata',      [R(1,1) screenD*R(1,1)/R(1,1)], 'ydata',[R(2,1) screenD*R(2,1)/R(1,1)], 'zdata',[R(3,1) screenD*R(3,1)/R(1,1)])
            set(app.Data.h.frame(2), 'xdata',      [n(1) screenD*n(1)/n(1)], 'ydata',[n(2) screenD*n(2)/n(1)], 'zdata',[n(3) screenD*n(3)/n(1)])




            set(app.Data.h.Vmeridian, 'xdata',-rinwVmeridianPoints(:,1), 'ydata', rinwVmeridianPoints(:,2), 'zdata', rinwVmeridianPoints(:,3));
            set(app.Data.h.Hmeridian, 'xdata',-rinwHmeridianPoints(:,1), 'ydata', rinwHmeridianPoints(:,2), 'zdata', rinwHmeridianPoints(:,3));

            set(app.Data.h.rVmeridianScreen, 'xdata',rinsVmeridianPoints(:,1), 'ydata', rinsVmeridianPoints(:,2), 'zdata', rinsVmeridianPoints(:,3) );
            set(app.Data.h.rHmeridianScreen, 'xdata',rinsHmeridianPoints(:,1), 'ydata', rinsHmeridianPoints(:,2), 'zdata', rinsHmeridianPoints(:,3) );
            set(app.Data.h.wVmeridianScreen, 'xdata',winwVmeridianPoints(:,1), 'ydata', winwVmeridianPoints(:,2), 'zdata', winwVmeridianPoints(:,3) );
            set(app.Data.h.wHmeridianScreen, 'xdata',winwHmeridianPoints(:,1), 'ydata', winwHmeridianPoints(:,2), 'zdata', winwHmeridianPoints(:,3) );

            
            [ListingPlaneX, ListingPlaneY,ListingPlaneZ] = Geometry3D.MakeListingsPlane(n);
            set(app.Data.h.listingsPlane, 'xdata',ListingPlaneX, 'ydata', ListingPlaneY, 'zdata',ListingPlaneZ);
            ax = Geometry3D.RotMat2AxisAngle(R);
            if Values.FollowListing_sLaw_ == "YES"
                ax2 = Geometry3D.RotMat2AxisAngle(RprimToV);
            else
                ax2 = [0 0 0]';
            end
            set(app.Data.h.rotationVector, 'xdata', -[0 ax(1)]*500, 'ydata', -[0 ax(2)]*500, 'zdata',-[0 ax(3)]*500);
            set(app.Data.h.rotationVectorFromPrimary, 'xdata', -[0 ax2(1)]*500, 'ydata', -[0 ax2(2)]*500, 'zdata',-[0 ax2(3)]*500);


            set(app.Data.h.image3Dplot.WorldFixedPath,  'xdata', wpatch3D(:,1), 'ydata',wpatch3D(:,2), 'zdata', wpatch3D(:,3));
            set(app.Data.h.image3Dplot.RetinaFixedPath,  'xdata', rpatch3D(:,1), 'ydata',rpatch3D(:,2), 'zdata', rpatch3D(:,3));



            % update screen plot

            set(app.Data.h.imageplot.wVmeridianScreen, 'xdata', winwVmeridianPoints(:,2)./winwVmeridianPoints(:,1)*screenD, 'ydata',winwVmeridianPoints(:,3)./winwVmeridianPoints(:,1)*screenD);
            set(app.Data.h.imageplot.wHmeridianScreen, 'xdata', winwHmeridianPoints(:,2)./winwHmeridianPoints(:,1)*screenD, 'ydata',winwHmeridianPoints(:,3)./winwHmeridianPoints(:,1)*screenD);
            set(app.Data.h.imageplot.rVmeridianScreen, 'xdata', rinwVmeridianPoints(:,2)./rinwVmeridianPoints(:,1)*screenD, 'ydata',rinwVmeridianPoints(:,3)./rinwVmeridianPoints(:,1)*screenD);
            set(app.Data.h.imageplot.rHmeridianScreen, 'xdata', rinwHmeridianPoints(:,2)./rinwHmeridianPoints(:,1)*screenD, 'ydata',rinwHmeridianPoints(:,3)./rinwHmeridianPoints(:,1)*screenD);
            
            Geometry3D.Update2DMeshGrid(app.Data.h.imageplot.screen, ygrid./xgrid*screenD, zgrid./xgrid*screenD)
            set(app.Data.h.imageplot.PrimaryPosition, 'xdata', n(2)./n(1)*screenD, 'ydata',n(3)./n(1)*screenD);
            set(app.Data.h.imageplot.FixationPosition,  'xdata', fy./fx*screenD, 'ydata',fz./fx*screenD);


            set(app.Data.h.imageplot.WorldFixedPath,  'xdata', wpatch(:,1), 'ydata',wpatch(:,2));
            set(app.Data.h.imageplot.RetinaFixedPath,  'xdata', rpatch(:,1), 'ydata',rpatch(:,2));
            % set(app.Data.h.imageplot.WorldFixedPathApprox,  'xdata', rinwpatchAprox(:,1), 'ydata',rinwpatchAprox(:,2));
            set(app.Data.h.imageplot.RetinaFixedPathApprox,  'xdata', rinwpatchAprox(:,1), 'ydata',rinwpatchAprox(:,2));


            % update retina plot

            set(app.Data.h.imageplot2.wVmeridianScreen, 'xdata', rad2deg(vmaz), 'ydata', rad2deg(-vmel));
            set(app.Data.h.imageplot2.wHmeridianScreen, 'xdata', rad2deg(-hmaz), 'ydata', rad2deg(hmel));
            set(app.Data.h.imageplot2.rVmeridianScreen, 'xdata',  [-0 0], 'ydata', [-90 90]);
            set(app.Data.h.imageplot2.rHmeridianScreen, 'xdata', [-90 90], 'ydata', [-0 0]);


            set(app.Data.h.imageplot2.WorldFixedPath,  'xdata', rad2deg(-wpatchaz), 'ydata', rad2deg(-wpatchel));
            set(app.Data.h.imageplot2.RetinaFixedPath,  'xdata',-rad2deg(rpatchaz), 'ydata', -rad2deg(rpatchel));
            set(app.Data.h.imageplot2.RetinaFixedPathApprox,  'xdata',-rad2deg(rpatchazApprox), 'ydata', -rad2deg(rpatchelApprox));
            %set(app.Data.h.imageplot2.WorldFixedPathApprox,  'xdata', rad2deg(winrpatchAprox(:,1)), 'ydata',rad2deg(winrpatchAprox(:,2)));

            meanError = mean(  sqrt( (rad2deg(rpatchaz) - rad2deg(rpatchazApprox)).^2 + (rad2deg(rpatchel) - rad2deg(rpatchelApprox)).^2));

            set(app.Data.h.imageplot2.ErrorText,'string',sprintf('Error from approximation = %0.1f deg \n', meanError))
        end
    end

    methods(Static) % Eye movement bascis

        function eyes = MakeEyes(ipdCm, fixationSpot, torsionVersion, torsionVergence, listingsPlanePitch, L2Gain, retinalShear, HHdeviation)
            
            if (~exist('listingsPlanePitch','var'))
                listingsPlanePitch = 0;
            end

            if (~exist('L2Gain','var'))
                L2Gain = 0.6;
            end

            if (~exist('retinalShear','var'))
                retinalShear = 0;
            end

            if (~exist('HHdeviation','var'))
                HHdeviation = 0;
            end

            % TODO: add troppias and fixation disparity

            % centers of the eyes
            eyes.R.Center = [0 -ipdCm/2 0]';
            eyes.L.Center = [0 ipdCm/2 0]';

            eyes.VergenceAngle = acosd(dot(normalize(fixationSpot - eyes.L.Center,'norm'), normalize(fixationSpot - eyes.R.Center,'norm')));

            % Primary position (orthogonal to listing's plane). Rotate the
            % plane according to the pitch first and the yaw according to
            % L2 
            % TODO: not sure if that relationship between pitch and yaw
            eyes.L.PrimaryPos = Geometry3D.RotZ(deg2rad(eyes.VergenceAngle) * L2Gain)  * Geometry3D.RotY(deg2rad(listingsPlanePitch)) * [1 0 0]';
            eyes.R.PrimaryPos = Geometry3D.RotZ(-deg2rad(eyes.VergenceAngle) * L2Gain) * Geometry3D.RotY(deg2rad(listingsPlanePitch)) * [1 0 0]';

            % Listing law rotations
            [listingsLeftRotMat, listingsLeftRotVec] = Geometry3D.LookAtListings(fixationSpot - eyes.L.Center, eyes.L.PrimaryPos);
            [listingsRightRotMat, listingsRightRotVec]  = Geometry3D.LookAtListings(fixationSpot - eyes.R.Center, eyes.R.PrimaryPos);

            % Torsion off of listing's plane
            torsionLeft = Geometry3D.RotX(deg2rad(torsionVersion - torsionVergence/2));
            torsionRight = Geometry3D.RotX(deg2rad(torsionVersion + torsionVergence/2));

            % Put everything together
            eyes.L.RotMat = listingsLeftRotMat * torsionLeft;
            eyes.R.RotMat = listingsRightRotMat * torsionRight;
            eyes.L.ListingsRotVec = listingsLeftRotVec;
            eyes.R.ListingsRotVec = listingsRightRotVec;

            % Retinal shear
            eyes.L.Shear = +retinalShear/2;
            eyes.R.Shear = -retinalShear/2;

            % Hering-Hilldebrand deviation
            eyes.L.HHdeviation = HHdeviation;
            eyes.R.HHdeviation = HHdeviation;

            % ------------

            % Build lines for drawing the eye

            % visual foveal directions
            lxfar = eyes.L.Center(1) + 10000*eyes.L.RotMat(1,1);
            lyfar = eyes.L.Center(2) + 10000*eyes.L.RotMat(2,1);
            lzfar = eyes.L.Center(3) + 10000*eyes.L.RotMat(3,1);
            rxfar = eyes.R.Center(1) + 10000*eyes.R.RotMat(1,1);
            ryfar = eyes.R.Center(2) + 10000*eyes.R.RotMat(2,1);
            rzfar = eyes.R.Center(3) + 10000*eyes.R.RotMat(3,1);

            % crosses around the eyes to illustrate torsion
            lrx = eyes.L.Center(1) - 10*eyes.L.RotMat(1,2);
            llx = eyes.L.Center(1) + 10*eyes.L.RotMat(1,2);
            lry = eyes.L.Center(2) - 10*eyes.L.RotMat(2,2);
            lly = eyes.L.Center(2) + 10*eyes.L.RotMat(2,2);
            lrz = eyes.L.Center(3) - 10*eyes.L.RotMat(3,2);
            llz = eyes.L.Center(3) + 10*eyes.L.RotMat(3,2);

            rrx = eyes.R.Center(1) - 10*eyes.R.RotMat(1,2);
            rlx = eyes.R.Center(1) + 10*eyes.R.RotMat(1,2);
            rry = eyes.R.Center(2) - 10*eyes.R.RotMat(2,2);
            rly = eyes.R.Center(2) + 10*eyes.R.RotMat(2,2);
            rrz = eyes.R.Center(3) - 10*eyes.R.RotMat(3,2);
            rlz = eyes.R.Center(3) + 10*eyes.R.RotMat(3,2);

            rbx = eyes.R.Center(1) - 10*eyes.R.RotMat(1,3);
            rtx = eyes.R.Center(1) + 10*eyes.R.RotMat(1,3);
            rby = eyes.R.Center(2) - 10*eyes.R.RotMat(2,3);
            rty = eyes.R.Center(2) + 10*eyes.R.RotMat(2,3);
            rbz = eyes.R.Center(3) - 10*eyes.R.RotMat(3,3);
            rtz = eyes.R.Center(3) + 10*eyes.R.RotMat(3,3);

            lbx = eyes.L.Center(1) - 10*eyes.L.RotMat(1,3);
            ltx = eyes.L.Center(1) + 10*eyes.L.RotMat(1,3);
            lby = eyes.L.Center(2) - 10*eyes.L.RotMat(2,3);
            lty = eyes.L.Center(2) + 10*eyes.L.RotMat(2,3);
            lbz = eyes.L.Center(3) - 10*eyes.L.RotMat(3,3);
            ltz = eyes.L.Center(3) + 10*eyes.L.RotMat(3,3);

            eyes.L.LinesForDrawing = ...
                [[eyes.L.Center(1) lxfar nan lbx ltx nan lrx llx] ; ...
                [eyes.L.Center(2) lyfar  nan lby lty nan lry lly] ; ...
                [eyes.L.Center(3) lzfar nan lbz ltz nan lrz llz] ]';
            eyes.R.LinesForDrawing = ...
                [[ eyes.R.Center(1) rxfar nan rbx rtx nan rrx rlx ] ; ...
                [ eyes.R.Center(2) ryfar nan rby rty nan rry rly ] ; ...
                [ eyes.R.Center(3) rzfar  nan rbz rtz nan rrz rlz] ]';

            [eyes.L.ListingPlaneX,eyes.L.ListingPlaneY,eyes.L.ListingPlaneZ] = Geometry3D.MakeListingsPlane(eyes.L.PrimaryPos);
            [eyes.R.ListingPlaneX,eyes.R.ListingPlaneY,eyes.R.ListingPlaneZ] = Geometry3D.MakeListingsPlane(eyes.R.PrimaryPos);
        end

        function [X,Y,Z] = MakeListingsPlane(n)

            % Normalize p to use it as the unit normal vector
            n = n / norm(n);

            % Choose a point on the plane (here, we use the origin)
            P0 = [0, 0, 0];

            % To define the plane, we need two orthonormal vectors (v1 and v2)
            % that lie in the plane (i.e., are perpendicular to n).

            % First, choose a vector 'a' that is not parallel to n
            if abs(n(1)) < abs(n(2)) && abs(n(1)) < abs(n(3))
                a = [1, 0, 0];
            elseif abs(n(2)) < abs(n(3))
                a = [0, 1, 0];
            else
                a = [0, 0, 1];
            end

            % Compute the first vector in the plane by taking the cross product
            v1 = cross(n, a);
            v1 = v1 / norm(v1);  % Normalize v1

            % Compute the second vector, which is perpendicular to both n and v1
            v2 = cross(n, v1);
            v2 = v2 / norm(v2);  % Normalize v2

            % Create a grid for the plane parameters (adjust the range as needed)
            [s, t] = meshgrid(-10:1:10, -10:1:10);

            % Parametric equation of the plane: P = P0 + s*v1 + t*v2
            X = P0(1) + s * v1(1) + t * v2(1);
            Y = P0(2) + s * v1(2) + t * v2(2);
            Z = P0(3) + s * v1(3) + t * v2(3);
        end

        function horopter = ComputeExtendedHoropter(eyes)

            % build the IDP vector pointing from the from the left to the
            % right eye
            ipdvector =  eyes.R.Center - eyes.L.Center;

            % First sample visual directions according to a helmholtz
            % coordinate system 2D grid
            [p, az, el] = Geometry3D.SampleVisualDirections(2000, 'Helmholtz', 90);

            % hering-hilldebrand deviation
            H = eyes.L.HHdeviation ;
            azlhh = 1/2 * (atan(2/H)  -acos(H./sqrt(H^2+4).*cos(2*az)));
            H = eyes.R.HHdeviation ;
            azrhh = 1/2 * (atan(2/H)  -acos(H./sqrt(H^2+4).*cos(2*az)));

            % shear the points according to the retinal shear
            azl = az + el*deg2rad(eyes.L.Shear) - azlhh;

            [x,y,z] = Geometry3D.HelmholtzToSphere(azl,el);
            pl = [x(:) y(:) z(:)];
            
            azr = az + el*deg2rad(eyes.R.Shear) + azrhh;
            
            [x,y,z] = Geometry3D.HelmholtzToSphere(azr,el);
            pr = [x(:) y(:) z(:)];

            % rotate the visual directions according to the eye position of
            % each eye
            pl = (eyes.L.RotMat*pl')';
            pr = (eyes.R.RotMat*pr')';

            H = zeros(size(p));
            D = zeros(size(p));
            ab = zeros(height(p),2);
            for i=1:height(pl)

                % solve the equations to find the coefficients of the
                % regression. ab tells you how for along the visual
                % direction of the left and the right eye the closest
                % corresponding point is. 
                ab(i,:) = pinv([pl(i,:)' -pr(i,:)'])*ipdvector;

                % calculate the 3D disparity of the points. That is the error
                % of the regression. 
                D(i,:) = [pl(i,:)' -pr(i,:)']*ab(i,:)' - ipdvector;

                % find the position of the point in space by averaging the
                % point of each eye
                hl = ab(i,1)*pl(i,:) + ipdvector'/2;
                hr = ab(i,2)*pr(i,:) - ipdvector'/2;
                H(i,:) = ( hl + hr  )/2;
            end

            % Calculate the length of the minimal disparity (in deg)
            D = 2*abs(atan2d( sqrt(sum(D.^2, 2)) , sqrt(sum(H.^2, 2))/2 ) );

            % Need to remove points where the crossing point is backwards
            D(ab(:,1)<0 | ab(:,2)<0,:) = nan;
            H(ab(:,1)<0 | ab(:,2)<0,:) = nan;

            % reshape to work with surf
            horopter = [];
            horopter.Points             = reshape(H, size(az,1), size(az,2), 3);
            horopter.MinDisparity       = reshape(D, size(az,1), size(az,2));
            horopter.VisualDirections   = p;
            horopter.Az                 = az;
            horopter.El                 = el;
        end

        function [R, rotvec, RfromPrimary] = LookAtListings(v, n)
            % [R, rotvec] = listingsLawRotation(v, n) 
            %
            % Computes the rotation matrix R that rotates the eye from
            % the reference position [1 0 0] to the direction v,
            % enforcing Listing's law according the primary position n.
            %
            % Another way of saying it. It gives the 3D pose of the eye
            % when looking at v by following Listing's law when the normal
            % to Listing's plane is n. 
            %            
            % Inputs:
            %   n : 3x1 vector, the normal to Listing's plane (primary position), (need not be normalized).
            %   v : 3x1 vector, the desired gaze direction (need not be normalized).
            %
            % Output:
            %   R : 3x3 rotation matrix of the pose of the eye looking at v
            %       this corresponds with a rotation from the identity
            %       matrix, not from the primary position.
            %   rotvec: rotation vector that corresponds with the rotation
            %       from primary position that is actually contained in Listing's plane
            %
            % The strategy is to compute the rotation from the primary
            % position to v, and from the primary position to the reference
            % position [1 0 0]. One corresponds with the actual rotation of
            % the eye and the other with a rotation of the coordinate
            % system.

            % 1. Normalize the target direction and primary positions
            v = normalize(v(:),'norm');  % ensure v is a column vector
            n = normalize(n(:),'norm');  % ensure n is a column vector

            % 2. Define the eye's reference position e1 = (1,0,0)
            r = [1; 0; 0];

            % 3. Compute the rotation axis using the cross product
            nv_ax  = cross(n, v);  
            nr_ax = cross(n, r); 

            % 4. Compute the rotation angle using the dot product
            nv_theta = acos(dot(n, v));
            nr_theta = acos(dot(n, r));

            % 5. Normalize the axis (check for degeneracy)
            nv_axnorm = norm(nv_ax);
            nr_axnorm = norm(nr_ax);
            
            %TODO: find a better way that does not need this check.
            nv_R = eye(3);
            nr_R = eye(3);

            if nv_axnorm > 1e-14
                nv_ax = nv_ax / nv_axnorm;

                % 7. Build the rotation matrix via Rodrigues' formula
                nv_R = Geometry3D.AxisAngle2Mat(nv_ax, nv_theta);
            end
            if nr_axnorm > 1e-14
                nr_ax = nr_ax / nr_axnorm;

                % 7. Build the rotation matrix via Rodrigues' formula
                nr_R = Geometry3D.AxisAngle2Mat(nr_ax, nr_theta);
            end

            % The final rotation is calculated as the rotation according to
            % listing's law on a reference frame aligned with the primary
            % position. And then rotated by the rotation of the primary
            % position. 
            R = nv_R*nr_R'; 
            
            rotvec = sin(nv_theta/2) * nv_ax;

            RfromPrimary =nv_R;
        end

        function [R] = LookAtListingsSimple(v)
            % listing's law rotation to look at v assuming that the primary
            % position is equal to the reference position [1 0 0]
            v = v/norm(v);
            x = v(1);
            y = v(2);
            z = v(3);
            d = 1+x;
            R = [x          -y                  -z; ...
                 y           1 - y.^2 ./ d      -y.*z ./ d; ...
                 z          -y.*z ./ d           1 - z.^2 ./ d ];
        end

        function [R, t] = LookAtCamera(eyeCenter, targetPos)
            % ---------------------------------------------------
            % Utility function: lookAtCamera
            % Creates an extrinsic [R,t] given camera center, target
            % ---------------------------------------------------

            upVec = [0 0 1]';
            forward = (eyeCenter - targetPos);
            forward = forward / norm(forward);

            right = cross(upVec, forward);
            right = right / norm(right);

            up = cross(forward, right);
            up = up / norm(up);

            R = [forward, right, up];  % 3x3 columns are forward, right, up, or row-based?

            % T is translation s.t. the camera center is at eyePos
            t =  R * eyeCenter(:);
        end

        function [eyepoints] = Points3DToEyes(points, eyes )

            ep = table();

            % Get the new xs in cm for the two eyes (transform the dots from head reference
            % to eye reference [nothing about eye angle here!])
            % Rotate the points to be in eye reference frame
            rPoints = (points{:,:}-eyes.R.Center')*eyes.R.RotMat;
            ep{:,{'RX' 'RY' 'RZ'  }} =  normalize(rPoints, 2, 'norm');
            lPoints = (points{:,:}-eyes.L.Center')*eyes.L.RotMat;
            ep{:,{'LX' 'LY' 'LZ'  }} =  normalize(lPoints, 2, 'norm');


            % shear the points according to the retinal shear
            [az, el] = Geometry3D.SphereToHelmholtz(ep.RX, ep.RY, ep.RZ);
            azl = az - el*deg2rad(eyes.R.Shear);
            [x,y,z] = Geometry3D.HelmholtzToSphere(azl,el);
            ep{:,{'RX' 'RY' 'RZ'  }} = [x(:) y(:) z(:)];
            [az, el] = Geometry3D.SphereToHelmholtz(ep.LX, ep.LY, ep.LZ);
            azr = az - el*deg2rad(eyes.L.Shear);
            [x,y,z] = Geometry3D.HelmholtzToSphere(azr,el);
            ep{:,{'LX' 'LY' 'LZ'  }} = [x(:) y(:) z(:)];




            % helmholtz coordinates make more sense for disparity
            % because the vertical rotation first does not change the
            % plane of regard
            [az, el] = Geometry3D.SphereToHelmholtz(ep.RX, ep.RY, ep.RZ);
            ep.RH = rad2deg(az);
            ep.RV = rad2deg(el);
            [az, el] = Geometry3D.SphereToHelmholtz(ep.LX, ep.LY, ep.LZ);
            ep.LH = rad2deg(az);
            ep.LV = rad2deg(el);

            % remove points past 90 deg in each direction
            badidx = abs(ep.RV)>=90 | abs(ep.RH)>=90 ;
            ep.RV(badidx) = nan;
            ep.RH(badidx) = nan;
            badidx = abs(ep.LV)>=90 | abs(ep.LH)>=90 ;
            ep.LV(badidx) = nan;
            ep.LH(badidx) = nan;

            % Get disparity
            ep.HDisparity = ep.LH - ep.RH;
            ep.VDisparity = ep.LV - ep.RV;

            eyepoints = ep;
        end

        function screen = MakeScreen(eyeCenter, distanceCm, sizeCm, resPix,  slantDeg)

            screen.distanceCm = distanceCm;
            screen.width = sizeCm(1);
            screen.height = sizeCm(2);
            screen.pixPerCmWidth = resPix(1)/sizeCm(1);
            screen.pixPerCmHeight = resPix(2)/sizeCm(2);
            screen.ScreenSlantDeg = slantDeg;



            % ---------------------------------------------------
            % 2) Define a camera via intrinsic & extrinsic
            % ---------------------------------------------------

            % 2a) Intrinsic matrix K
            %   
            fx = resPix(1)/sizeCm(1)*screen.distanceCm*cosd(slantDeg);
            fy = resPix(2)/sizeCm(2)*screen.distanceCm*cosd(slantDeg);
            cx = resPix(1)/2 - sind(slantDeg)*resPix(1)/sizeCm(1);
            cy = resPix(2)/2;
            K = [cx, fx,   0;
                cy,  0,   -fy;
                1    0,   0 ];

            % TODO: slant does not work

            % 2b) Define extrinsic from a known camera "eye" position & orientation
            %
            % Let's place the camera at EYE_POS in the world, and aim it along -Z
            % with Y up, etc. The standard "lookAt" approach can help define R.

            % We'll compute R,t by a simple "lookAt" utility 
            screenNormalToEye = [distanceCm*cosd(slantDeg)*cosd(slantDeg);distanceCm*sind(slantDeg)*cosd(slantDeg);0];
            [R, t] = Geometry3D.LookAtCamera(eyeCenter, eyeCenter + screenNormalToEye);

            % rotate according to the slant
            %R = Geometry3D.RotZ(deg2rad(slantDeg))*R;

            % Now the extrinsic matrix is [R | t].

            % Full camera projection matrix
            screen.ProjectionMatrix = K * [R, t];
        end

        function screenPoints  = PointsEyesToScreen(worldPoints, screen)

            % ---------------------------------------------------
            % 3) Project 3D world points
            % ---------------------------------------------------

            % Convert 3D worldPoints to homogeneous coords: [X Y Z 1]
            homog3D = [worldPoints{:,:}, ones(height(worldPoints),1)]';  % 4xN

            % Project: x = P * X_world
            % Then do perspective division
            proj = screen.ProjectionMatrix * homog3D;   % 3xN
            uv   = proj ./ proj(3,:);  % each column: [u; v; 1]

            % Transpose for Nx2
            screenPoints = uv(1:2, :)';  % Nx2
        end
    end

    methods(Static) % BASIC 3D coordinate conversion functions (all in radians)

        % Single axis rotation matrices

        function R = RotZ(theta)
            % RotZ: Returns the rotation matrix about the Z-axis. 3rd axis
            % in a right handed coordinate system.
            %   Input: theta - rotation angle in radians.
            %   (Haslwanter 1995, eq. 2)

            R = [ cos(theta)  -sin(theta)   0;
                sin(theta)   cos(theta)   0;
                0            0            1];
        end

        function R = RotY(phi)
            % RotY: Returns the rotation matrix about the Y-axis. 2nd axis
            % in a right handed coordinate system.
            %   Input: phi - rotation angle in radians.
            %   (Haslwanter 1995, eq. 3)

            R = [ cos(phi)   0    sin(phi);
                0          1    0;
                -sin(phi)   0    cos(phi)];
        end

        function R = RotX(psi)
            % RotX: Returns the rotation matrix about the X-axis. 1st axis
            % in a right handed coordinate system.
            %   Input: psi - rotation angle in radians.
            %   (Haslwanter 1995, eq. 4)

            R = [ 1     0          0;
                0     cos(psi)  -sin(psi);
                0     sin(psi)   cos(psi)];
        end

        % Euler angles to rotation matrix conversions

        function R = EulerAngles2RotMat(HVT,coordinateSystem)
            switch(coordinateSystem)
                case 'Fick'
                    R = Geometry3D.Fick2RotMat(HVT);
                case 'Helmholtz'
                    R = Geometry3D.Helm2RotMat(HVT);
                case 'Harms'
                    R = Geometry3D.Harms2RotMat(HVT);
                case 'Hess'
                    R = Geometry3D.Hess2RotMat(HVT);
                case 'ImagePlane'
                    R = Geometry3D.ImagePlane2RotMat(HVT);
                case 'Polar'
                    R = Geometry3D.Polar2RotMat(HVT);
            end
        end

        function HVT = RotMat2EulerAngles(R,coordinateSystem)
            switch(coordinateSystem)
                case 'Fick'
                    HVT = Geometry3D.RotMat2Fick(R);
                case 'Helmholtz'
                    HVT = Geometry3D.RotMat2Helm(R);
                case 'Harms'
                    HVT = Geometry3D.RotMat2Harms(R);
                case 'Hess'
                    HVT = Geometry3D.RotMat2Hess(R);
                case 'ImagePlane'
                    HVT = Geometry3D.RotMat2ImagePlane(R);
                case 'Polar'
                    HVT = [0 0 0]; % TODO
            end
        end

        function R = Fick2RotMat(HVT)

            H = HVT(1);
            V = HVT(2);
            T = HVT(3);

            R = Geometry3D.RotZ(H) * Geometry3D.RotY(V) * Geometry3D.RotX(T);
        end
        
        function HVT = RotMat2Fick(M)

            r31 = M(3,1);
            r21 = M(2,1);
            r32 = M(3,2);

            HVT(2) = -asin(r31);
            HVT(1) = asin(r21/cos(HVT(2)));
            HVT(3) = asin(r32/cos(HVT(2)));
        end

        function R = Helm2RotMat(HVT)

            H = HVT(1);
            V = HVT(2);
            T = HVT(3);

            R = Geometry3D.RotY(V) * Geometry3D.RotZ(H) * Geometry3D.RotX(T);
        end

        function HVT = RotMat2Helm(M)
            r21 = M(2,1);
            r31 = M(3,1);
            r23 = M(2,3);

            HVT(1) = asin(r21);
            HVT(2) = -asin(r31/cos(HVT(1)));
            HVT(3) = -asin(r23/cos(HVT(1)));
        end

        function R = Hess2RotMat(HVT)
            % For this coordinate system torsion is not well defined
            % we will use rotation with axis perpendicular to x and then
            % add torsin

            [x, y, z] = Geometry3D.HessToSphere(HVT(1), -HVT(2));
            R = Geometry3D.LookAtListingsSimple([x y z])*Geometry3D.RotX(HVT(3))';
        end
        function HVT = RotMat2Hess(M)
            % For this coordinate system torsion is not well defined
            % we will use rotation with axis perpendicular to x and then
            % add torsin
            [az, el] = Geometry3D.SphereToHess(M(:,1));

            R = Geometry3D.LookAtListingsSimple(M(:,1))*M';

            T = acos(R(1,1));

            HVT(1) = az;
            HVT(2) = el;
            HVT(3) = T;
        end

        function R = Harms2RotMat(HVT)
            % For this coordinate system torsion is not well defined
            % we will use rotation with axis perpendicular to x and then
            % add torsin

            [x, y, z] = Geometry3D.HarmsToSphere(HVT(1), -HVT(2));
            R = Geometry3D.LookAtListingsSimple([x y z])*Geometry3D.RotX(HVT(3))';
        end
        function HVT = RotMat2Harms(M)
            % For this coordinate system torsion is not well defined
            % we will use rotation with axis perpendicular to x and then
            % add torsin

            [az, el] = Geometry3D.SphereToHarms(M(:,1));

            R = Geometry3D.LookAtListingsSimple(M(:,1))*M';

            T = acos(R(1,1));

            HVT(1) = az;
            HVT(2) = el;
            HVT(3) = T;
        end

        function R = ImagePlane2RotMat(HVT,f)
            if (~exist('f','var'))
                f = 1;
            end
            % For this coordinate system torsion is not well defined
            % we will use rotation with axis perpendicular to x and then
            % add torsin

            [x, y, z] = Geometry3D.ImagePlaneToSphere(HVT(1), -HVT(2),f);
            R = Geometry3D.LookAtListingsSimple([x y z])*Geometry3D.RotX(HVT(3))';
        end
        function HVT = RotMat2ImagePlane(M,f)
            if (~exist('f','var'))
                f = 1;
            end
            % For this coordinate system torsion is not well defined
            % we will use rotation with axis perpendicular to x and then
            % add torsin

            [az, el] = Geometry3D.SphereToImagePlane(M(1,1), M(2,1), M(3,1),f);

            R = Geometry3D.LookAtListingsSimple(M(:,1))*M';

            T = acos(R(1,1));

            HVT(1) = az;
            HVT(2) = el;
            HVT(3) = T;
        end

        function R = Polar2RotMat(AET)

            % the input is angle, eccentricity and torsion
            % sequency of rotations is as follows
            % 1- rotate the coordinate system according the angle (
            % then rotate an eccentricity ammount arond the rotated Z axis
            % finally undo the intial torsion and add the actual torsion
            A = AET(1);
            E = AET(2);
            T = AET(3);

            R = Geometry3D.RotX(A) * Geometry3D.RotZ(E) *Geometry3D.RotX(-A)*Geometry3D.RotX(T);
        end

        % Quaternion and rotation vectors to rotation matrix conversions

        function q = RotMat2Quat(M)
            % From quaternion navy book
            M = M';
            m11 = M(1,1);
            m12 = M(1,2);
            m21 = M(2,1);
            m22 = M(2,2);
            m33 = M(3,3);
            m23 = M(2,3);
            m32 = M(3,2);
            m31 = M(3,1);
            m13 = M(1,3);

            q0 = 1/2*sqrt(1+m11+m22+m33);

            q = [q0 (m23-m32)/(4*q0) (m31-m13)/(4*q0) (m12-m21)/(4*q0)];
        end

        function R = Quat2RotMat(q)
            % from quaternion dynamics pdf

            E = [...
                -q(2) q(1) -q(4) q(3); ...
                -q(3) q(4) q(1) -q(2); ...
                -q(4) -q(3) q(2) q(1); ...
                ];

            G = [...
                -q(2) q(1) q(4) -q(3); ...
                -q(3) -q(4) q(1) q(2); ...
                -q(4) q(3) -q(2) q(1); ...
                ];

            R = E*G';

        end

        function [axis, angle] = Quat2AxisAngle(q)
            angle = 2*acos(q(1));
            if ( angle ~= 0)
                axis = (q(2:4)/sin(angle));
                axis = axis/norm(axis);
            else
                axis = [1 0 0];
            end
        end

        function q = AxisAngle2Quat(axis, angle)
            q = [cos(angle/2) sin(angle/2)*axis/norm(axis)];
        end

        function R = AxisAngle2Mat(axis, angle)
            % Helper function to build a rotation matrix
            % from a unit rotation axis and a rotation angle.
            x = axis(1);
            y = axis(2);
            z = axis(3);

            c = cos(angle);
            s = sin(angle);
            C = 1 - c;

            % Rodrigues' rotation formula
            R = [ ...
                c + x*x*C,   x*y*C - z*s,  x*z*C + y*s;
                y*x*C + z*s, c + y*y*C,    y*z*C - x*s;
                z*x*C - y*s, z*y*C + x*s,  c + z*z*C   ];
        end

        function [axis, angle] = RotMat2AxisAngle(R)
            % HASLWANTER (1995) equation 23
            r = (1/ (1+R(1,1)+R(2,2)+R(3,3)) ) * ([R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)]);

            norm_r = norm(r);
            axis = r/norm_r;
            angle = asin(norm_r)/2;
        end


        % Functions to converte between spherical coordinates and coordinate
        % systems for 2D rotations. When going from 2D to 3D it always
        % gives normalized vectors on the unit sphere.

        function [az, el] = SphereToAzEl(x,y,z, coordinateSystem)
            switch(coordinateSystem)
                case 'Fick'
                    [az, el] = Geometry3D.SphereToFick(x,y,z);
                case 'Helmholtz'
                    [az, el] = Geometry3D.SphereToHelmholtz(x,y,z);
                case 'Harms'
                    [az, el] = Geometry3D.SphereToHarms(x,y,z);
                case 'Hess'
                    [az, el] = Geometry3D.SphereToHess(x,y,z);
                case 'ImagePlane'
                    [az, el] = Geometry3D.SphereToImagePlane(x,y,z);
                case 'Polar'
                    [az, el] = Geometry3D.SphereToPolar(x,y,z);
            end
        end

        function [x,y,z] = AzElToSphere(az,el, coordinateSystem)
            switch(coordinateSystem)
                case 'Fick'
                    [x,y,z] = Geometry3D.FickToSphere(az, el);
                case 'Helmholtz'
                    [x,y,z] = Geometry3D.HelmholtzToSphere(az, el);
                case 'Harms'
                    [x,y,z] = Geometry3D.HarmsToSphere(az, el);
                case 'Hess'
                    [x,y,z] = Geometry3D.HessToSphere(az, el);
                case 'ImagePlane'
                    [x,y,z] = Geometry3D.ImagePlaneToSphere(az, el);
                case 'Polar'
                    [x,y,z] = Geometry3D.PolarToSphere(az, el);
            end
        end


        function [x, y, z] = FickToSphere(az, el)
            x = cos( el ) .* cos( az );
            y = sin( az ) .* cos( el );   % longitude
            z = sin( el );                % latitude
        end

        function [az,el] = SphereToFick(x,y,z)
            D = sqrt( x.^2 + y.^2 + z.^2 );

            az = atan2( y, x );   % longitude
            el = asin( z ./ D );   % latitude
        end

        function [x, y, z] = HelmholtzToSphere(az,el)

            x = cos( az ) .* cos( el );
            y = sin( az );                 % latitude
            z = sin( el ) .* cos( az );   % longitude
        end

        function [az,el] = SphereToHelmholtz(x,y,z)
            D = sqrt( x.^2 + y.^2 + z.^2 );

            az = asin( y ./ D );  % latitude
            el = atan2( z, x );  % longitude
        end

        function [x, y, z] = HessToSphere(az, el)

            z = sin(el);
            y = sin(az);
            x = sqrt( 1 - z.^2 - y.^2 );

            % Need to flip negative azimuts
            x(cos(az)<0) = -x(cos(az)<0);

            % Need to force zero at 90 deg azimuth to avoid some complex
            % numbers that can appear numerically
            x(z.^2 + y.^2 == 1) = 0;

            % some of the combinations of azimuth and elevation actually
            % don't exist within the sphere.
            outsidepoints = (z.^2 + y.^2) > 1;
            z(outsidepoints) = nan;
            y(outsidepoints) = nan;
            x(outsidepoints) = nan;
        end

        function [az,el] = SphereToHess(x,y,z)
            D = sqrt( x.^2 + y.^2 + z.^2 );

            az = asin( y ./ D ); % latitude
            el = asin( z ./ D ); % latitude
        end

        function [x, y, z] = HarmsToSphere(az, el)
            x = 1 ./ sqrt(1 + tan(az).^2 + tan(el).^2 );

            
            x(cos(az)==0 | cos(el)==0) = 0;

            y = tan(az) .* x;
            z = tan(el) .* x;


            % z(abs(cos(el))<eps & sin(el) < 0 ) = sin(az(abs(cos(el))<eps  & sin(el)  < 0));
            % y(abs(cos(el))<eps & sin(el)  > 0 ) = cos(az(abs(cos(el))<eps  & sin(el)  > 0));
            % 
            % z(abs(cos(el))<eps & sin(el)  > 0 ) = sin(az(abs(cos(el))<eps  & sin(el)  > 0));
            % y(abs(cos(el))<eps & sin(el)  < 0 ) = -cos(az(abs(cos(el))<eps  & sin(el)  < 0));
            % 
            % 
            % z(abs(cos(az))<eps & sin(az) < 0 ) = sin(el(abs(cos(az))<eps  & sin(az)  < 0));
            % y(abs(cos(az))<eps & sin(az)  > 0 ) = -cos(el(abs(cos(az))<eps  & sin(az)  > 0));
            % 
            % z(abs(cos(az))<eps & sin(az)  > 0 ) = sin(el(abs(cos(az))<eps  & sin(az)  > 0));
            % y(abs(cos(az))<eps & sin(az)  < 0 ) = cos(el(abs(cos(az))<eps  & sin(az)  < 0));
        end

        function [az,el] = SphereToHarms(x,y,z)
            az = atan2(y,x);  % longitude
            el = atan2(z,x);  % longitude
        end

        function [x,y,z] = ImagePlaneToSphere(u,v, f)
            if (~exist('f','var'))
                f = 1;
            end
            denom = sqrt(f.^2 + u.^2 + v.^2);
            x = f ./ denom;
            y = u ./ denom;
            z = v ./ denom;
        end

        function [u, v] = SphereToImagePlane(x,y,z,f)
            if (~exist('f','var'))
                f = 1;
            end
            u = f*y./x;
            v = f*z./x;
        end

        function [x, y, z] = PolarToSphere(angle, eccentricity)
            % the input is angle, eccentricity 
            % sequency of rotations is as follows
            % 1- rotate the coordinate system around the x axis according the angle
            % then rotate an eccentricity ammount arond the rotated Z axis
            x = cos( eccentricity );
            y = sin( eccentricity ) .* cos( angle );
            z = sin( eccentricity ) .* sin( angle );
        end

        function [angle, ecc] = SphereToPolar(x,y,z)
            ecc = acos(x);
            angle = atan2(z,y);
        end

    end

    methods(Static) % Geometric projections

        function projPoints = StereographicProjectionZY(spherePoints)
            % StereographicProjectionZY computes the stereographic projection of points
            % on a sphere onto the zy-plane (x = 1) using the projection point (-1,0,0).
            %
            % Input:
            %   spherePoints - an N-by-3 matrix where each row is a point [x, y, z]
            %                  that lies on the sphere centered at the origin.
            %
            % Output:
            %   projPoints   - an N-by-3 matrix of projected points [0, y', z'] where:
            %                  y' = y/(1-x) and z' = z/(1-x)
            %
            % Note:
            %   Points with x = 1 (or very close to 1) will produce a division by zero
            %   and are set to NaN.

            % Extract the coordinates
            x = spherePoints(:, 1);
            y = spherePoints(:, 2);
            z = spherePoints(:, 3);

            % Compute the denominator (1-x)
            denom =  x+1;

            % Check for points where denom is nearly zero
            if any(abs(denom) < 1e-10)
                warning('Some points have x ~ 1, mapping to infinity; their projection will be NaN.');
            end

            % Compute the projected y and z coordinates
            y_proj = 2*y ./ denom;
            z_proj = 2*z ./ denom;

            % The projected x coordinate is 0 for the zy-plane
            x_proj = zeros(size(x));

            % Combine into the output matrix
            projPoints = [x_proj, y_proj, z_proj];
        end

        function projPoints = AzimuthalEquidistantProjectionZY(spherePoints)
            % AzimuthalEquidistantProjectionZY projects points on the unit sphere onto the
            % tangent plane at the north pole (x = 1) using the azimuthal equidistant projection.
            %
            % Input:
            %   spherePoints - an N-by-3 matrix where each row is a point [x, y, z] on the unit sphere.
            %
            % Output:
            %   projPoints   - an N-by-3 matrix of projected points [X, Y, 1] on the plane z = 1.
            %
            % The projection is computed as follows:
            %   r     = acos(x)         % Angular distance from the north pole
            %   theta = atan2(y, z)       % Azimuth angle in the xy-plane
            %   X     = r * cos(theta)
            %   Y     = r * sin(theta)
            %
            % Note:
            %   This projection maps the north pole (0,0,1) to (0,0,1) and is defined for all points
            %   on the sphere.

            % Extract coordinates from the input
            x = spherePoints(:,1);
            y = spherePoints(:,2);
            z = spherePoints(:,3);

            % Compute the angular distance from the north pole (in radians)
            r = acos(x);

            % Compute the azimuth angle
            theta = atan2(z, y);

            % Compute the projected coordinates in the tangent plane at the north pole
            Y_proj = r .* cos(theta);
            Z_proj = r .* sin(theta);

            % The projected x-coordinate is 1 (on the plane x = 1)
            projPoints = [ones(size(x)), Y_proj, Z_proj];
        end
    end

    methods(Static) % Graphic helpers
        function [f, ax] = setup3DPlot(axlim)
            f = figure('color','w');
            hold on

            set(gca, 'visible','off')
            axis([-axlim axlim -axlim axlim -axlim axlim])
            xlim([-axlim axlim])
            ylim([-axlim axlim])
            zlim([-axlim axlim])


            line([-axlim axlim], [0 0], [0 0],'color',[0.6 0.6 0.6])
            line([0 0], [-axlim axlim], [0 0],'color',[0.6 0.6 0.6])
            line([0 0], [0 0], [-axlim axlim],'color',[0.6 0.6 0.6])
            text(axlim*1.2,0,0,'x')
            text(0,axlim*1.2,0,'y')
            text(0,0,axlim*1.2,'z')
            axis square
            view(-45,25)
            ax = gca;
        end

        function displayRefFrame(Xnew,Ynew,Znew)
            set(gca,'nextplot','add')
            if ( nargin == 1)
                M = Xnew;
                Xnew = M(:,1);
                Ynew = M(:,2);
                Znew = M(:,3);
            end
            quiver3(0,0,0,Xnew(1),Xnew(2),Xnew(3), 'linewidth',2)
            quiver3(0,0,0,Ynew(1),Ynew(2),Ynew(3), 'linewidth',2)
            quiver3(0,0,0,Znew(1),Znew(2),Znew(3), 'linewidth',2)
            text(Xnew(1),Xnew(2),Xnew(3),'r_x')
            text(Ynew(1),Ynew(2),Ynew(3),'r_y')
            text(Znew(1),Znew(2),Znew(3),'r_z')
        end

        function DrawPolarGrid(R, numCircles, numRadialLines, labelAngle, addlabels )
            if ( ~exist('addlabels','var'))
                addlabels = 1;
            end
            theta = linspace(0, 2*pi, 360);     % Angle values for circles

            % Plot concentric circles
            radii = linspace(R/numCircles, R, numCircles);
            for r = radii
                xCircle = r * cos(theta);
                yCircle = r * sin(theta);
                plot(xCircle, yCircle, 'k:');  % Use a dotted line style
            end
            % Add text labels for each circle indicating the distance from the center
            for r = radii
                xLabel = r * cos(labelAngle);
                yLabel = r * sin(labelAngle);
                text(xLabel, yLabel, sprintf('%.0f', r), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
            end
            if ( addlabels)
                text(90, 0, {'World' 'right'}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
                text(-90, 0, {'World' 'left'}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
                text(0, 90, 'World up', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
                text(0, -90, 'World down', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            end

            % Plot radial lines
            angles = linspace(0, 2*pi, numRadialLines+1); % +1 to complete the circle
            for a = angles
                xLine = [0, R * cos(a)];
                yLine = [0, R * sin(a)];
                plot(xLine, yLine, 'k:');      % Use a dotted line style
            end
        end

        function h = Plot2DMeshGrid(X,Y)
            X(end+1,:) = nan;
            Y(end+1,:) = nan;
            X(:,end+1) = nan;
            Y(:, end+1) = nan;
            XX = [X;X'];
            YY = [Y;Y'];
            h = plot(XX(:),YY(:));
        end

        function Update2DMeshGrid(h,X,Y)
            X(end+1,:) = nan;
            Y(end+1 ,:)= nan;
            X(:,end+1) = nan;
            Y(:, end+1) = nan;
            XX = [X;X'];
            YY = [Y;Y'];
            set(h,'xdata',XX(:),'ydata',YY(:));
        end
    end

    %% Spherical coordinate methods
    methods(Static)
        function PlotTangentSphereJacobians()
            %%
            %
            R = 1; % radius of the eye
            step = 1;
            range = 90;
            [az, el] = meshgrid(-2*range:step:2*range,-range:step:range); % azimuths and elevations to include

            [xs,ys,zs] = Geometry3D.FickToSphere(deg2rad(az),deg2rad(el));
            a = ParticleSampleSphere( 'N' ,200);
            x = a(a(:,1)>0.01,1);
            y = a(a(:,1)>0.01,2);
            z = a(a(:,1)>0.01,3);
            v = [x(:), y(:), z(:) ];


            [Jacobian.rotational.a, Jacobian.linear.a, ~, Jacobian3D.rotational.a, Jacobian3D.linear.a] = Geometry3D.CalculateMotionJacobianFields(v,'TangentSphere');
            duwxdt = -squeeze(Jacobian.rotational.a(1,1,:));
            duwydt = -squeeze(Jacobian.rotational.a(1,2,:));
            duwzdt = -squeeze(Jacobian.rotational.a(1,3,:));
            dvwxdt = -squeeze(Jacobian.rotational.a(2,1,:));
            dvwydt = -squeeze(Jacobian.rotational.a(2,2,:));
            dvwzdt = -squeeze(Jacobian.rotational.a(2,3,:));

            dudx = -squeeze(Jacobian.linear.a(1,1,:));
            dudy = -squeeze(Jacobian.linear.a(1,2,:));
            dudz = -squeeze(Jacobian.linear.a(1,3,:));
            dvdx = -squeeze(Jacobian.linear.a(2,1,:));
            dvdy = -squeeze(Jacobian.linear.a(2,2,:));
            dvdz = -squeeze(Jacobian.linear.a(2,3,:));

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

            fac = 10;
            hl = [];
            hl(1) = quiver3(x,y,z, dudx(:)/fac, dudy(:)/fac, dudz(:)/fac,'linewidth',2,'DisplayName', '\partial xyz / \partial u','AutoScale','off')
            hl(2) = quiver3(x,y,z, dvdx(:)/fac, dvdy(:)/fac, dvdz(:)/fac,'linewidth',2,'DisplayName', '\partial xyz / \partial v','AutoScale','off')
            % quiver3(x,y,z, x, y, z,'linewidth',1)

            legend(hl,'Location', 'northeast','box','off','fontsize',14);




            % figure('color','w')
            subplot(2,5,3,'nextplot','add')
            axis equal;
            set(gca,'xlim',[-1 1], 'ylim',[-1 1])
            quiver(y.*acosd(x)/90,z.*acosd(x)/90, dudx, dvdx,'linewidth',2)
            xlabel( '\partial u / \partial x','fontsize',14);
            ylabel( '\partial v / \partial x','fontsize',14);
            set(gca,'xtick',[],'ytick',[])

            subplot(2,5,4,'nextplot','add')
            axis equal;
            set(gca,'xlim',[-1 1], 'ylim',[-1 1])
            quiver(y.*acosd(x)/90,z.*acosd(x)/90, dudy, dvdy,'linewidth',2)
            xlabel( '\partial u / \partial y','fontsize',14);
            ylabel( '\partial v / \partial y','fontsize',14);
            set(gca,'xtick',[],'ytick',[])

            subplot(2,5,5,'nextplot','add')
            axis equal;
            set(gca,'xlim',[-1 1], 'ylim',[-1 1])
            quiver(y.*acosd(x)/90,z.*acosd(x)/90, dudz, dvdz,'linewidth',2)
            xlabel( '\partial u / \partial z','fontsize',14);
            ylabel( '\partial v / \partial z','fontsize',14);
            set(gca,'xtick',[],'ytick',[])


            subplot(2,5,8,'nextplot','add')
            axis equal;
            set(gca,'xlim',[-1 1], 'ylim',[-1 1])
            quiver(y.*acosd(x)/90,z.*acosd(x)/90, duwxdt, dvwxdt,'linewidth',2)
            xlabel( '\partial u / \omega_x dt','fontsize',14);
            ylabel( '\partial v / \omega_x dt','fontsize',14);
            set(gca,'xtick',[],'ytick',[])

            subplot(2,5,9,'nextplot','add')
            axis equal;
            set(gca,'xlim',[-1 1], 'ylim',[-1 1])
            quiver(y.*acosd(x)/90,z.*acosd(x)/90, duwydt, dvwydt,'linewidth',2)
            xlabel( '\partial u / \omega_y dt','fontsize',14);
            ylabel( '\partial v / \omega_y dt','fontsize',14);
            set(gca,'xtick',[],'ytick',[])

            subplot(2,5,10,'nextplot','add')
            axis equal;
            set(gca,'xlim',[-1 1], 'ylim',[-1 1])
            quiver(y.*acosd(x)/90,z.*acosd(x)/90, duwzdt, dvwzdt,'linewidth',2)
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
            [x,y,z] = Geometry3D.FickToSphere(deg2rad([0;50]),deg2rad([0;40]));

            v = [x(:), y(:), z(:) ];
            [Jacobian.rotational.a, Jacobian.linear.a, ~, Jacobian3D.rotational.a, Jacobian3D.linear.a] = Geometry3D.CalculateMotionJacobianFields(v,'TangentSphere');

            dudx = -squeeze(Jacobian.linear.a(1,1,:));
            dudy = -squeeze(Jacobian.linear.a(1,2,:));
            dudz = -squeeze(Jacobian.linear.a(1,3,:));
            dvdx = -squeeze(Jacobian.linear.a(2,1,:));
            dvdy = -squeeze(Jacobian.linear.a(2,2,:));
            dvdz = -squeeze(Jacobian.linear.a(2,3,:));


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
            hl(1) = quiver3(x,y,z, dudx, dudy, dudz,'linewidth',2,'DisplayName', 'u');
            hl(2) = quiver3(x,y,z, dvdx, dvdy, dvdz,'linewidth',2,'DisplayName', 'v');

            legend(hl,'Location', 'northeast','box','off','fontsize',14);
        end

        function PlotMotionTemplatesForSphericalCoordinateSystems(type)
            %% Visualization of different spherical coordinate systems and motion flows projected onto them
            %
            % x y z is a right handed reference frame
            %
            % The eye points towards x and up is z
            %
            %

            % clear all, close all

            if ( ~exist('type','var'))
                type = 'analytical';
            end

            step = 10;
            range = 80;
            [az, el] = meshgrid(-range:step:range,-range:step:range); % azimuths and elevations to include
            coordinateSystems = {'Fick', 'Helmholtz', 'Harms','Hess'};
            colorsflow = {'r' 'b' 'k'};
            flows = {'linear' 'rotational'};
            dims = {'x' 'y' 'z'};

            % vectors and axis
            m.linear.x = [1 0 0 0.3 0 0];
            m.linear.y = [0 1 0 0 0.3 0];
            m.linear.z = [0 0 1 0 0 .3];
            m.rotational.x = [-2 0 0 2 0 0];
            m.rotational.y = [0 -2 0 0 2 0];
            m.rotational.z = [0 0 -2 0 0 2];

            f = figure('color','white');
            tiledlayout(5,7,"TileSpacing","tight","Padding","tight");


            f.Position = [f.Position(1) f.Position(2) 4*f.Position(3) f.Position(4)];

            % Do the 3D flow plots on the top row
            v = Geometry3D.SampleVisualDirections(numel(az));
            [~, ~, ~, Jacobian3D.rotational.a, Jacobian3D.linear.a] = Geometry3D.CalculateMotionJacobianFields(v);

            nexttile

            for iFlow = 1:2
                J3D = Jacobian3D.(flows{iFlow}).(lower(type(1)));

                for vi=1:3

                    ax = nexttile;
                    title(['Flow due to ' flows{iFlow}  ' motion along ',dims{vi}])
                    set(ax,'nextplot','add')
                    view(140,15)
                    axis equal;
                    [xs,ys,zs] = sphere(100);
                    mesh(xs,ys,zs,'FaceAlpha', 0.5,'facecolor',0.8*[1 1 1],'EdgeColor','none');

                    quiver3(v(:,1),v(:,2),v(:,3),squeeze(rad2deg(J3D(1,vi,:))),squeeze(rad2deg(J3D(2,vi,:))),squeeze(rad2deg(J3D(3,vi,:))), 'linewidth',1.5)

                    if ( iFlow==1)
                        quiver3(m.linear.(dims{vi})(1), m.linear.(dims{vi})(2), m.linear.(dims{vi})(3), m.linear.(dims{vi})(4), m.linear.(dims{vi})(5), m.linear.(dims{vi})(6) , 'color','r', 'linewidth',2)
                    else
                        line(m.linear.(dims{vi})([1 4]), m.linear.(dims{vi})([2 5]), m.linear.(dims{vi})([3 6]), 'color','r','linewidth',2,'linestyle','-.')
                    end
                    set(gca,'xlim',[-1.3 1.3],'ylim',[-1.3 1.3],'zlim',[-1.3 1.3])

                end
            end

            % Draw coordinate systems and flows for each system in each row

            for i=1:length(coordinateSystems)
                sys = coordinateSystems{i};

                nexttile
                set(gca, 'nextplot','add')

                % calculate spherical coordinates depending on the
                % coordinate system to make the grids
                switch(lower(sys))
                    case 'fick'
                        [x,y,z] = Geometry3D.FickToSphere(deg2rad(az),deg2rad(el));
                    case 'helmholtz'
                        [x,y,z]  = Geometry3D.HelmholtzToSphere(deg2rad(az),deg2rad(el));
                    case 'harms'
                        [x,y,z] = Geometry3D.HarmsToSphere(deg2rad(az),deg2rad(el));
                    case 'hess'
                        [x,y,z] = Geometry3D.HessToSphere(deg2rad(az),deg2rad(el));
                end
                
                % Calculate jacobians
                v = [x(:) y(:) z(:)];
                [Jacobian.rotational.a, Jacobian.linear.a, ~, Jacobian3D.rotational.a, Jacobian3D.linear.a] = Geometry3D.CalculateMotionJacobianFields(v,sys);
                [Jacobian.rotational.n, Jacobian.linear.n, ~, Jacobian3D.rotational.n, Jacobian3D.linear.n] = Geometry3D.CalculateMotionJacobianFieldsNumerical(v,sys);

                % check if numerical and analytical jacobians match
                if ( mean((Jacobian.rotational.a - Jacobian.rotational.n).^2 , "all", 'omitnan') / mean((Jacobian.rotational.a).^2 +(Jacobian.rotational.n).^2 , "all", 'omitnan') > 1 )
                    disp(['Numerical and analytical jacobians do not match for ' flows{iFlow}  ' motion along ',dims{vi}])
                end
                if ( mean((Jacobian.linear.a - Jacobian.linear.n).^2,"all" , 'omitnan') / mean((Jacobian.linear.a).^2 +(Jacobian.linear.n).^2 , "all", 'omitnan') > 1 )
                    disp(['Numerical and analytical jacobians do not match for ' flows{iFlow}  ' motion along ',dims{vi}])
                end

                % draw sphere in first column
                view(125,15)
                axis equal;
                mesh(x,y,z ,'FaceAlpha', 0.9,'facecolor',[1 1 1]);
                xlabel(gca,'x'), ylabel(gca,'y'), zlabel(gca,'z')
                % line of sight
                line([0 1.5 ],[0 0 ],[0 0 ],'color','r','linewidth',2)
                title(sys)

                % draw template jacobians in other columns
                for iFlow = 1:2
                    J = Jacobian.(flows{iFlow}).(lower(type(1)));

                    for vi=1:3
                        nexttile
                        set(gca,'nextplot','add')
                        grid on

                        quiver(az(:), el(:), squeeze(rad2deg(J(1,vi,:))), squeeze(rad2deg(J(2,vi,:))), colorsflow{vi}, 'linewidth',1.5)

                        plot(0,0,'ro')


                        axis equal;
                        set(gca,'xlim',[-92 92],'ylim',[-92 92])
                        set(gca,'xtick',[-90:30:90],'ytick',[-90:30:90])
                    end
                end
            end
            delete(nexttile(1))
        end

        function J = FickInverseJacobian(az, el)

            % dx/daz, dy/daz, dz/daz
            J(1,1,:) = -sin( az ) .* cos( el );
            J(2,1,:) = cos( az ) .* cos( el );
            J(3,1,:) = 0;

            % dx/del, dy/del, dz/del
            J(1,2,:) = -cos( az ) .* sin( el );
            J(2,2,:) = -sin( az ) .* sin( el );
            J(3,2,:) = cos( el );
        end

        function J  = HelmholtzInverseJacobian(az, el)

            % dx/daz, dy/daz, dz/daz
            J(1,1,:) = -sin(az) .* cos(el);
            J(2,1,:) = cos(az);
            J(3,1,:) = -sin(el) .* sin(az);

            % dx/del, dy/del, dz/del
            J(1,2,:) = -cos(az) .* sin(el);
            J(2,2,:) = zeros(size(x));
            J(3,2,:) = cos(az) .* cos(el);
        end

        function J = HarmsInverseJacobian(az, el)
            % Define the trigonometric functions
            tan_az = tan(az);
            sec_az = sec(az);
            tan_el = tan(el);
            sec_el = sec(el);

            % Common denominator
            denominator = (1 + tan_az.^2 + tan_el.^2).^(3/2);

            % Calculate each element of the Jacobian matrix

            % dx/daz, dy/daz, dz/daz
            J(1,1,:) = -tan_az .* sec_az.^2 ./ denominator;
            J(2,1,:) = (sec_az.^2 - tan_az.^2 .* sec_az.^2) ./ denominator;
            J(3,1,:) = -tan_az .* sec_az.^2 .* tan_el ./ denominator;
            
            % dx/del, dy/del, dz/del
            J(1,2,:) = -tan_el .* sec_el.^2 ./ denominator;
            J(2,2,:) = -tan_az .* sec_az.^2 .* tan_el ./ denominator;

            J(3,2,:) = (sec_el.^2 - tan_el.^2 .* sec_el.^2) ./ denominator;

        end

        function [x, y, z] = SpiralSphere(Ni)

            N=round(Ni(1)*2.1);

            gr=(1+sqrt(5))/2;       % golden ratio
            ga=2*pi*(1-1/gr);       % golden angle

            i=(0:(N-1))';              % particle (i.e., point sample) index
            lat=acos(1-2*i/(N-1));  % latitude is defined so that particle index is proportional to surface area between 0 and lat
            lon=i*ga;               % position particles at even intervals along longitude

            % Convert from spherical to Cartesian coordinates
            x=sin(lat).*cos(lon);
            y=sin(lat).*sin(lon);
            z=cos(lat);

            % ensure N points in the positive side of the sphere only
            [~,idx] = sort(x);
            idx = idx(end-(Ni-1):end);
            x = x(idx);
            y = y(idx);
            z = z(idx);
        end

        function [x,y,z] = RandomSampleSphere(Ni)

            % Partition the [-1,1]x[0,2*pi] domain into ceil(sqrt(N))^2 subdomains
            % and then draw a random sample for each
            n=ceil(sqrt(Ni));
            ds=2/n;
            [Xc,Yc]=meshgrid((-1+ds/2):ds:(1-ds/2));

            x=ds*(rand(n^2,1)-0.5);
            y=ds*(rand(n^2,1)-0.5);

            x=x+Xc(:);
            y=y+Yc(:);
            clear Xc Yc

            % Remove excess samples
            R=n^2-Ni;
            if R>0
                idx=randperm(n^2,R);
                x(idx)=[];
                y(idx)=[];
            end

            lon=(x+1)*pi;
            z=y;

            % Convert z to latitude
            z(z<-1)=-1;
            z(z>1)=1;
            lat=acos(z);

            % Convert spherical to rectangular co-ords
            s=sin(lat);
            x=cos(lon).*s;
            y=sin(lon).*s;

        end

    end

    methods(Static) % misc

        function PlotLinearCalibrationError()
            trueAngle = -50:50;
            calibrationAngle = 15;
            estimatedAngle = calibrationAngle/(tand(calibrationAngle))*tand(trueAngle);

            figure
            subplot(1,2,1)
            plot(trueAngle,trueAngle);
            hold
            plot(trueAngle,estimatedAngle)

            plot(calibrationAngle*[-1 0 1], calibrationAngle*[-1 0 1],'o')

            xlabel('True angle')
            ylabel('Estimated angle with linear calibration.')

            subplot(1,2,2)
            plot(trueAngle,trueAngle-trueAngle);
            hold
            plot(trueAngle,estimatedAngle-trueAngle)

            plot(calibrationAngle*[-1 0 1], 0*calibrationAngle*[-1 0 1],'o')

            xlabel('True angle')
            ylabel('Estimated angle error with linear calibration.')

            set(gca,'ylim',[-5 5])

        end
    end
end
