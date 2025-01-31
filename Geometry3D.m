classdef Geometry3D
    %Geometry3D Summary of this class goes here
    %   Detailed explanation goes here
    %
    %   Coordinate system
    %       For 3D world
    %       Z is positive forward
    %       X is positive to the right
    %       system)
    %       Y is positive up
    %
    %       For eye reference frame
    %       Z is positive up
    %       X is positive forward
    %       system)
    %       Y is positive left
    %
    %       For euler angles we use helmholtz
    %       H is positive towards the right
    %       V is positive upwards
    %       T is positive top to the right
    %
    %       Disparity is L - R
    %

    properties
    end

    methods(Static)


        function demoCoordinateSystemsAndPlane()

            app = InteractiveUI('Coordinate systems',@(app) (Geometry3D.demoCoordinateSystemsAndPlaneUpdate(app)), 0.1);
            app.AddDropDown('Coordinate System',   1,  [ "TangentSphere", "Fick", "Helmholtz", "Harms","Hess"])
            % app.AddSlider('Eye Radius',           0.02,  [0.01    1])
            app.AddSlider('Azimuth',           -20,  [-90 90])
            app.AddSlider('Elevation',           20,  [-90 90])
            % app.AddSlider('Torsion Version',      0,  [-20  20])
            % app.AddSlider('Torsion Vergence',     0,  [-20  20])
            %app.AddSlider('Ground plane slant',          0,  [-90  90])
            % app.AddSlider('Ground plane tilt',           0,  [0    90])
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

        function visualDirections = SampleVisualDirections(N, type)
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


            if ( ~exist('type','var'))
                type = 'Spiral';
            end

            % get the sample visual directions in spherical coordinates
            % depending on the coordinate system
            N = round(sqrt(N)).^2; % make sure N is a square number
            range = 80;
            step = range*2/(sqrt(N)-1);
            [az, el] = meshgrid(deg2rad(-range:step:range),deg2rad(-range:step:range)); % azimuths and elevations to include

            switch(type)
                case 'Fick'
                    [x,y,z] = Geometry3D.FickToSphere(az,el);
                case 'Helmholtz'
                    [x,y,z] = Geometry3D.HelmholtzToSphere(az,el);
                case 'Harms'
                    [x,y,z] = Geometry3D.HarmsToSphere(az,el);
                case 'Hess'
                    [x,y,z] = Geometry3D.HessToSphere(az,el);
                case 'Spiral'
                    [x, y, z] = Geometry3D.SpiralSphere(N);
                case 'Random'
                    [x, y, z] = Geometry3D.RandomSampleSphere(N);
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
            [Jv, Jw] = CalculateMotionJacobianFields([0,0,0]);
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

        function [Jw, Jv] = CalculateMotionJacobianFields(visualDirections, coordSys)
            if ( ~exist('coordSys','var'))
                coordSys = 'TangentSphere';
            end

            N = height(visualDirections);

            switch(coordSys)
                case 'Fick'
                    [az, el] = Geometry3D.SphereToFick(visualDirections(:,1),visualDirections(:,2),visualDirections(:,3));
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.FickLinearJacobian(az, el);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.FickRotationalJacobian(az, el);
                case 'Helmholtz'
                    [az, el] = Geometry3D.SphereToHelmholtz(visualDirections(:,1),visualDirections(:,2),visualDirections(:,3));
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HelmholtzLinearJacobian(az, el);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.HelmholtzRotationalJacobian(az, el);
                case 'Harms'
                    [az, el] = Geometry3D.SphereToHarms(visualDirections(:,1),visualDirections(:,2),visualDirections(:,3));
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HarmsLinearJacobian(az, el);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.HarmsRotationalJacobian(az, el);
                case 'Hess'
                    [az, el] = Geometry3D.SphereToHess(visualDirections(:,1),visualDirections(:,2),visualDirections(:,3));
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HessLinearJacobian(az, el);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.HessRotationalJacobian(az, el);
                case 'TangentSphere'
                    x = visualDirections(:,1);
                    y = visualDirections(:,2);
                    z = visualDirections(:,3);
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.TangentSphereLinearJacobian(x,y,z);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.TangentSphereRotationalJacobian(x,y,z);
            end

            % collect the jacobians into 3D matrices 2x3xN
            % (2x3 jacobian at each visual direciton)
            Jv = cat(1, reshape([dazdx(:)'; dazdy(:)'; dazdz(:)'],1,3, N),  reshape([deldx(:)'; deldy(:)' ;deldz(:)'],1,3, N));
            Jw = cat(1, reshape([rdazwx(:)'; rdazwy(:)'; rdazwz(:)'],1,3, N),  reshape([rdelwx(:)' ;rdelwy(:)' ;rdelwz(:)'],1,3, N));
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

        function hs = DisplayMotionField(motionField, visualDirections, hs)


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
            visualDirections = Geometry3D.SampleVisualDirections(N,coordSys);
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
            switch(coordSys)
                case 'Fick'
                    [px,py,pz] = Geometry3D.FickToSphere(paz,pel);
                    [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.FickLinearInverseJacobian(paz, pel);
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.FickLinearJacobian(paz, pel);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.FickRotationalJacobian(paz, pel);

                case 'Helmholtz'
                    [px,py,pz] = Geometry3D.HelmholtzToSphere(paz,pel);
                    [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.HelmholtzLinearInverseJacobian(paz, pel);
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HelmholtzLinearJacobian(paz, pel);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.HelmholtzRotationalJacobian(paz, pel);
                case 'Harms'
                    [px,py,pz] = Geometry3D.HarmsToSphere(paz,pel);
                    [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.HarmsLinearInverseJacobian(paz, pel);
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HarmsLinearJacobian(paz, pel);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.HarmsRotationalJacobian(paz, pel);
                case 'Hess'
                    [px,py,pz] = Geometry3D.HessToSphere(paz,pel);
                    [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.HessLinearInverseJacobian(paz, pel);
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HessLinearJacobian(paz, pel);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.HessRotationalJacobian(paz, pel);
                case 'TangentSphere'
                    [px,py,pz] = Geometry3D.FickToSphere(paz,pel);
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.TangentSphereLinearJacobian(px,py,pz);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.TangentSphereRotationalJacobian(px,py,pz);
                    [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.TangentSphereInverseJacobian(px,py,pz);
            end

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

        function demoCoordinateSystems(WHICHFLOW)

            if ( nargin <1)
                Geometry3D.demoCoordinateSystems('linear')
                Geometry3D.demoCoordinateSystems('rotational')
                return
            end

            R = 0.025; % radius of the eye
            step = 10;
            [az, el] = meshgrid(-80:step:80,-80:step:80); % azimuths and elevations to include
            coordinateSystems = { 'Listings', 'Fick','Helmholtz', 'Harms','Hess'};




            f = figure('color','white');
            f.Position = [f.Position(1) f.Position(2) 4*f.Position(3) f.Position(4)];

            for i=1:length(coordinateSystems)
                sys = coordinateSystems{i};
                subplot(4,length(coordinateSystems),i,'nextplot','add');

                % calculate spherical coordinates depending on the coordinate system
                switch(sys)
                    case 'Fick'
                        [x,y,z] = Geometry3D.FickToSphere(az,el);
                    case 'Helmholtz'
                        [x,y,z] = Geometry3D.HelmholtzToSphere(az,el);
                    case 'Listings'
                        [x,y,z] = Geometry3D.ListingsToSphere(az,el);
                    case 'Harms'
                        [x,y,z] = Geometry3D.HarmsToSphere(az,el);
                    case 'Hess'
                        [x,y,z] = Geometry3D.HessToSphere(az,el);
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
                % This this placements actually needs to be multiplied by
                % the ration between the eye radius and the distance of the
                % object seen by this receptor. Right now wea assume ration
                % of 1.
                dw = 0.1; % differential angular rotation in deg

                az2 = az;
                el2 = el;

                colorsflow = {'r' 'b' 'k'};
                dims = {'x' 'y' 'z'};

                dxyz = {[dv,0,0],[0,dv,0],[0,0,dv]};
                dRM = {Geometry3D.RotX(dw), Geometry3D.RotY(dw), Geometry3D.RotZ(dw)};

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
                            [az2, el2 ] = Geometry3D.SphereToFick(x2,y2,z2);
                        case 'Helmholtz'
                            [az2, el2 ] = Geometry3D.SphereToHelmholtz(x2,y2,z2);
                        case 'Listings'
                            az2 = nan(size(az));
                            el2 = nan(size(el));
                        case 'Harms'
                            [az2, el2 ] = Geometry3D.SphereToHarms(x2,y2,z2);
                        case 'Hess'
                            [az2, el2 ] = Geometry3D.SphereToHess(x2,y2,z2);
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
        end

        function demoDisparity()

            app = InteractiveUI('Disparity Simulator',@(app) (Geometry3D.demoDisparityUpdate(app)), 0.2);
            app.AddDropDown('Stimulus',         1,  ["CROSS" "RANDOMPLANE" "GRIDPLANE" "HLINE" "VLINE"])
            app.AddSlider('IPD mm',             60, [10 100])
            app.AddSlider('Stimulus Scale',     1,  [0.1 10])
            app.AddSlider('Stimulus Distance',  40, [10 200])
            app.AddSlider('Fixation Distance',  30, [10 200])
            app.AddSlider('Fixation X',         0,  [-100 100])
            app.AddSlider('Fixation Y',         0,  [-100 100])
            app.AddSlider('Torsion Version',    0,  [-20 20])
            app.AddSlider('Torsion Vergence',   0,  [-20 20])
            app.AddSlider('Plane slant',        0,  [-90 90])
            app.AddSlider('Plane Tilt',         0,  [0 90])
            app.AddDropDown('View3D',           1,  ["Oblique" "TOP" "SIDE"])
            app.AddSlider('Screen slant',       0,  [-30 30])
            app.AddDropDown('Simulate torsion', 1,  ["Yes" "No"])


            app.Data.Screen = struct();
            app.Data.Screen.SizeCm = [30*16/9 30];
            app.Data.Screen.ResPix = [1920 1080];
            app.Data.Screen.Distance = 57;
            app.Data.Screen.Slant = 0;

            app.Data.FixationSpot = struct();
            app.Data.FixationSpot.X = 0;
            app.Data.FixationSpot.Y = 0;
            app.Data.FixationSpot.Z = 0;

            app.Open();

        end

        function exmapleCalculateEyePoints()
            %%

            SizeCm = [30*16/9 30];
            ResPix = [1920 1080];
            ScreenDistance = 57;
            ScreenSlant = 0;

            IPDMm = 60;
            FixationSpot = [];
            FixationSpot.X = 0;
            FixationSpot.Y = 0;
            FixationSpot.Z = 50; % fixation spot distance
            TorsionVersion = 0;
            TorsionVergence = 0;

            StimulusDistance = 50;
            sizeStimCm = 40;
            PlaneTilt = 0;
            PlaneSlant = 0;

            % make world points
            numDots = 200;
            worldPoints = table();
            worldPoints.X = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
            worldPoints.Y = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
            worldPoints.Z = zeros(numDots, 1);

            worldPoints.X = worldPoints.X;
            worldPoints.Y = worldPoints.Y;
            worldPoints.Z = zeros(size(worldPoints.Z));
            % Rotate by slant and tilt
            R = Geometry3D.Quat2RotMat(Geometry3D.AxisAngle2Quat([cosd(PlaneTilt) sind(PlaneTilt) 0],deg2rad(PlaneSlant)));
            worldPoints{:,:} = (R*worldPoints{:,:}')';
            % Displace by distance
            worldPoints.Z = worldPoints.Z + StimulusDistance;
            % end make world points



            leftEyeScreen = Geometry3D.MakeScreen(SizeCm, ResPix, ScreenDistance, ScreenSlant);
            rightEyeScreen = Geometry3D.MakeScreen(SizeCm, ResPix, ScreenDistance, ScreenSlant);
            eyes = Geometry3D.MakeEyes(IPDMm/10, FixationSpot, TorsionVersion, TorsionVergence);

            eyePoints = Geometry3D.Points3DToEyes(worldPoints, eyes);
            screenPoints = Geometry3D.PointsEyesToScreen(eyes, eyePoints, leftEyeScreen, rightEyeScreen);
        end

        function demoMotionFlow()

            app = InteractiveUI('Motion Flow Simulator',@(app) (Geometry3D.demoMotionFlowUpdate(app)), 0.2);
            app.AddDropDown('Stimulus',      1,  ["Ground plane" "Point cloud"])
            app.AddSlider('Eye height',           100,[0    500])
            % app.AddSlider('Stimulus Scale',       1,  [0.1  10])
            app.AddSlider('Fixation Distance',    30, [10   200])
            app.AddSlider('Fixation X',           0,  [-100 100])
            app.AddSlider('Fixation Y',           0,  [-100 100])
            % app.AddSlider('Torsion Version',      0,  [-20  20])
            % app.AddSlider('Torsion Vergence',     0,  [-20  20])
            app.AddSlider('Ground plane slant',          0,  [-90  90])
            % app.AddSlider('Ground plane tilt',           0,  [0    90])
            app.AddSlider('Angular velocity X',   0,  [-100 100])
            app.AddSlider('Angular velocity Y',   0,  [-100 100])
            app.AddSlider('Angular velocity Z',   0,  [-100 100])
            app.AddSlider('Linear velocity X',    0,  [-100 100])
            app.AddSlider('Linear velocity Y',    0,  [-100 100])
            app.AddSlider('Linear velocity Z',    0,  [-100 100])
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

        function demoQuaternion()
            app = InteractiveUI('Quaternion demo',@(app) (Geometry3D.demoQuaternionUpdate(app)), 1);
            app.AddSlider('q0',1, [-1 1])
            app.AddSlider('q1',0, [-1 1])
            app.AddSlider('q2',0, [-1 1])
            app.AddSlider('q3',0, [-1 1])

            % app.AddDropDown('2D HVT Coordinate system',      1,  ["Helmholtz" "Fick" "Harns" "Hess"])
            app.AddDropDown('2D HVT Coordinate system',      1,  ["Helmholtz" "Fick" "Harms" "Hess"])
            app.AddSlider('H',0, [-90 90])
            app.AddSlider('V',0, [-90 90])
            app.AddSlider('T',0, [-90 90])

            % app.AddSlider('Listings angle',0, [-90 90])
            % app.AddSlider('Listings eccentricity',0, [-90 90])
            % app.AddSlider('Listings torsion',0, [-90 90])

            % app.AddSlider('x',0, [-1 1])
            % app.AddSlider('y',0, [-1 1])
            % app.AddSlider('z',0, [-1 1])
            % app.AddSlider('a',0, [-180 180])


            app.Open();
        end

        function eyes = MakeEyes(ipdCm, fixationSpot, torsionVersion, torsionVergence)

            % TODO: add deviations. Not sure how yet. Probably need six
            % numbers to add all posible deviations. It would probably also
            % change assumptions about the fixation spot and how it is
            % being drawn

            % Position of eyes
            eyes.R.X = ipdCm/2;
            eyes.R.Y = 0;
            eyes.R.Z = 0;
            eyes.L.X = - ipdCm/2;%
            eyes.L.Y = 0;
            eyes.L.Z = 0;

            % Orientation of eyes in helmoltz coordinates, they make more
            % sense for disparity because it naturally encodes first the
            % plane of regard. That is, the plane containing the fixation
            % spot and the center of the eyes. The horizontal eye position
            % is the angle within that plane.
            %
            % right is positive, Up is positive, top-pole to right is
            % positive
            eyes.R.V = atan2d(fixationSpot.Y-eyes.R.Y, fixationSpot.Z);
            eyes.R.H = atan2d(fixationSpot.X-eyes.R.X, fixationSpot.Z./(cosd(eyes.R.V)));
            eyes.R.T = torsionVersion - torsionVergence;

            eyes.L.V = atan2d(fixationSpot.Y-eyes.L.Y, fixationSpot.Z);
            eyes.L.H = atan2d(fixationSpot.X-eyes.L.X, fixationSpot.Z./(cosd(eyes.L.V)));
            eyes.L.T = torsionVersion + torsionVergence;


            % Get the rotation matrices describing the eye angles in
            % helmholtz order
            % Right hand rotations around Y are to the left, so we flip the
            % sign of the rotation
            RH = Geometry3D.RotZ(-deg2rad(eyes.R.H));
            RV = Geometry3D.RotY(deg2rad(eyes.R.V));
            RT = Geometry3D.RotX(deg2rad(eyes.R.T));
            eyes.R.RM = RV*RH*RT; % Rot matrix defining the orientation of the eye relative to straight ahead direction.

            LH = Geometry3D.RotZ(-deg2rad(eyes.L.H));
            LV = Geometry3D.RotY(deg2rad(eyes.L.V));
            LT = Geometry3D.RotX(deg2rad(eyes.L.T));
            eyes.L.RM = LV*LH*LT; % Rot matrix defining the orientation of the eye relative to straight ahead direction.

        end

        function screen = MakeScreen(sizeCm, resPix, distance, slant)

            screen.widthCm = sizeCm(1);
            screen.heightCm = sizeCm(2);
            screen.pixPerCmWidth = resPix(1)/sizeCm(1);
            screen.pixPerCmHeight = resPix(2)/sizeCm(2);
            screen.middleX = resPix(1)/2;
            screen.middleY = resPix(2)/2;

            screen.ScreenDistance = distance;
            screen.ScreenDistance = distance;

            screen.ScreenSlant = slant;
        end


        function [eyepoints] = Points3DToEyes(points, eyes )


            ep = table();

            % The eye coordinate system is right handed
            % with x pointing out
            % y pointing to left
            % z pointing up
            %
            % it is a bit messy to have a different coordinate system than
            % for the 3D scene, but I can't think of the rotations any
            % other way
            %

            % Get the new xs in cm for the two eyes (transform the dots from head reference
            % to eye reference [nothing about eye angle here!])
            ep.RY0 = -(points.X - eyes.R.X);
            ep.RZ0 = -points.Y;
            ep.RX0 = points.Z;
            ep.LY0 = -(points.X - eyes.L.X);
            ep.LZ0 = -points.Y;
            ep.LX0 = points.Z;

            % Rotate the points to be in eye reference frame
            ep{:,{'RX' 'RY' 'RZ'  }} =  ep{:,{'RX0' 'RY0' 'RZ0'  }}*eyes.R.RM;
            ep{:,{'LX' 'LY' 'LZ'  }} =  ep{:,{'LX0' 'LY0' 'LZ0'  }}*eyes.L.RM;

            % helmholtz coordinates make more sense for disparity
            % because the vertical rotation first does not change the
            % plane of regard
            ep.RV = atan2d(ep.RZ, ep.RX);
            ep.RH = atan2d(-ep.RY, ep.RX./abs(cosd(ep.RV)));
            ep.LV = atan2d(ep.LZ, ep.LX);
            ep.LH = atan2d(-ep.LY, ep.LX./abs(cosd(ep.LV)));

            % remove points past 90 deg in each direction
            badidx = abs(ep.RV)>=90 | abs(ep.RH)>=90 ;
            ep.RV(badidx) = nan;
            ep.RH(badidx) = nan;
            badidx = abs(ep.LV)>=90 | abs(ep.LH)>=90 ;
            ep.LV(badidx) = nan;
            ep.LH(badidx) = nan;

            % Get disparity!! % TODO, mayber this should be 3D disparity
            ep.HDisparity = ep.LH - ep.RH;
            ep.VDisparity = ep.LV - ep.RV;

            eyepoints = ep;
        end


        function screenPoints  = PointsEyesToScreen(eyes, eyePoints, LeyeScreen, ReyeScreen)

            % TODO: to simulate torsion we want to rotate the points only
            % around the axis that goes trhough the fixation spot by how
            % much the torsion is.
            % This should work off of the points shifted (the ones with
            % zero at the end). That way we don't need to undo the
            % horizontal and vertical rotation that was done to rotate the
            % points to full eye coordinates

            screenPoints.LX = LeyeScreen.middleX + LeyeScreen.ScreenDistance*tand(eyes.L.H + eyePoints.LH)*LeyeScreen.pixPerCmWidth;
            screenPoints.LY = LeyeScreen.middleY + LeyeScreen.ScreenDistance*tand(-eyes.L.V + eyePoints.LV)*LeyeScreen.pixPerCmHeight;
            screenPoints.RX = ReyeScreen.middleX + ReyeScreen.ScreenDistance*tand(eyes.R.H + eyePoints.RH)*ReyeScreen.pixPerCmWidth;
            screenPoints.RY = ReyeScreen.middleY + ReyeScreen.ScreenDistance*tand(-eyes.R.V + eyePoints.RV)*ReyeScreen.pixPerCmHeight;
        end


        function demo()
            %%
            Rx = Geometry3D.RotX(deg2rad(0));
            Ry = Geometry3D.RotY(deg2rad(0));
            Rz = Geometry3D.RotZ(deg2rad(30));

            x  = [1 0 0];
            y  = [0 1 0];
            z  = [0 0 1];

            rx = x*Rx*Ry*Rz;
            ry = y*Rx*Ry*Rz;
            rz = z*Rx*Ry*Rz;


            axlim = 1.5;
            [f, ax] = Geometry3D.setup3DPlot(axlim);
            Geometry3D.displayRefFrame(rx,ry,rz)
            %%
        end
    end

    methods(Static, Access = private) % DEMO DISPARITY

        function [f, heyes, hfix, hscreen, hpoints, hspoints, hLRpoints, hdisparity, hLscreen] = demoDisparityInitPlots(points, eyePoints, screenPoints, eyes, leftEyeScreen, rightScreen)

            % View the dots, eyes, and fixation dot
            scr_siz = get(0,'ScreenSize');
            margin = floor(0.1*(scr_siz(4)));
            f = figure('color','w','position',floor([margin margin scr_siz(3)*2.8/4 scr_siz(4)*2/4 ]));
            axes('OuterPosition',[ 0    0    0.5    1], 'nextplot','add')

            hpoints = line(0,0, 0,'linestyle','none','marker','o','Color','k');
            hspoints = line(0,0, 0,'linestyle','none','marker','o','Color',[0.8 0.8 0.8]);
            hfix = line(0,0, 0,'linestyle','none','marker','o','Color','r','LineWidth',2, 'markersize',20);

            heyes.c(1) = plot3([eyes.L.X], [eyes.L.Y ], 0,'o','Color','b', 'markersize',10); % left eye fixation spot and right eye
            heyes.c(2) = plot3([eyes.R.X], [eyes.R.Y], 0,'o','Color','r', 'markersize',10); % left eye fixation spot and right eye
            heyes.l = line(0,0,0,'linestyle','-','Color','b'); % left eye fixation spot and right eye
            heyes.r = line(0,0,0,'linestyle','-','Color','r'); % left eye fixation spot and right eye

            grid
            xlim([-100 100]), zlim([-100 100]), ylim([-5 200]);
            xlabel('X (cm)'), zlabel('Y (cm)'), ylabel('Z (cm)');
            view(-60,10)
            title('3D world');
            hscreen = [];
            % hscreen(1) = plot3([-1 -1 1 1 -1]*leftEyeScreen.widthCm/2+eyes.L.X, [1 1 1 1 1]*leftEyeScreen.ScreenDistance, [-1 1 1 -1 -1]*leftEyeScreen.heightCm/2);
            % hscreen(2) = plot3([-1 -1 1 1 -1]*rightScreen.widthCm/2+eyes.R.X, [1 1 1 1 1]*rightScreen.ScreenDistance, [-1 1 1 -1 -1]*rightScreen.heightCm/2);
            hLRpoints = [];
            t = eyePoints;

            polaraxes('OuterPosition',[0.44    0.3    0.28    0.7], 'nextplot','add')
            hLRpoints.L = polarplot(zeros(size(t.LV)), zeros(size(t.LV)),'o');
            hLRpoints.R = polarplot(zeros(size(t.RV)), zeros(size(t.RV)),'o');
            hLRpoints.FP = polarplot(0,0,'ro','linewidth',3);
            grid
            % xlim([min([t.RH t.LH],[],"all") max([t.RH t.LH],[],"all")])
            % ylim([min([t.RV t.LV],[],"all") max([t.RV t.LV],[],"all")])
            % xlabel('X (deg)')
            % ylabel('Y (deg)')
            set(gca,"RLim",[0 90])
            set(gca,'ThetaTickLabel',[])
            grid
            % legend({'Right Eye','Left Eye'})
            title('Visual field Left eye (blue) Right eye (red)');

            %-----------------------------
            % disparity plot
            %-----------------------------


            axes('OuterPosition',[ 0.65    0.3    0.35    0.7], 'nextplot','add')
            set(gca, 'PlotBoxAspectRatio',[1 1 1])
            hdisparity = quiver(t.RH, t.RV, t.LH-t.RH, t.LV-t.RV, 'AutoScale', "off");
            grid
            title('Disparity (Helmholtz, L->R)')
            set(gca,'xlim',[-40 40],'ylim',[-40 40])

            % subplot(2,4,7)
            % for i = 1:height(t)
            %     line([t.RH(i,1) t.LH(i,1)],[t.RV(i,1) t.LV(i,1)])
            % end
            % title('Disparity With Torsion')


            axes('OuterPosition',[0.5493    0.0091    0.1566    0.3412], 'nextplot','add')
            hLscreen.LPoints = plot(screenPoints.LX(1:end-1), screenPoints.LY(1:end-1), '+');

            hLscreen.LFP = plot(screenPoints.LX(end), screenPoints.LY(end), 'ro');
            grid
            set(gca,'xlim',[0 1920],'ylim',[0 1080])
            set(gca,'PlotBoxAspectRatio',[16 9 1])
            title('Left eye screen (sim torsion)');



            axes('OuterPosition',[  0.7625    0.0165    0.1566    0.3412], 'nextplot','add')
            hLscreen.RPoints = plot(screenPoints.RX(1:end-1), screenPoints.RY(1:end-1), '+');

            hLscreen.RFP = plot(screenPoints.RX(end), screenPoints.RY(end), 'ro');
            grid
            set(gca,'xlim',[0 1920],'ylim',[0 1080])
            set(gca,'PlotBoxAspectRatio',[16 9 1])
            title('Right eye screen (sim torsion)');

        end

        function demoDisparityUpdate(app)

            Values = app.Values;

            if (~isfield(app.Data,"stimulus"))
                app.Data.stimulus = "NONE";
            end

            %%
            sizeStimCm = 40;

            if (app.Data.stimulus ~= string(Values.Stimulus) )
                % Make some dots in world coordinates (XYZ)
                %
                % Make a table to start putting everything together
                switch (Values.Stimulus)
                    case 'RANDOMPLANE'
                        numDots = 200;
                        worldPoints = table();
                        worldPoints.X = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
                        worldPoints.Y = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
                        worldPoints.Z = zeros(numDots, 1);
                    case 'GRIDPLANE'
                        % TODO: need to update the update funciton so the Z is not overwritten by a
                        % constant distance. Just shift things.
                        [X,Y] = meshgrid(-(sizeStimCm/2):2:(sizeStimCm/2),-(sizeStimCm/2):2:(sizeStimCm/2));
                        worldPoints = table();
                        worldPoints.X = X(:);
                        worldPoints.Y = Y(:);
                        worldPoints.Z = zeros(size(worldPoints.X));

                    case 'HLINE'
                        % TODO: need to update the update funciton so the Z is not overwritten by a
                        % constant distance. Just shift things.
                        [X,Y] = meshgrid(-(sizeStimCm/2):0.5:(sizeStimCm/2),(0)*ones(size(-(sizeStimCm/2):0.5:(sizeStimCm/2))));
                        worldPoints = table();
                        worldPoints.X = X(:);
                        worldPoints.Y = Y(:);
                        worldPoints.Z = zeros(size(worldPoints.X));
                    case 'VLINE'
                        % TODO: need to update the update funciton so the Z is not overwritten by a
                        % constant distance. Just shift things.
                        [X,Y] = meshgrid((0)*ones(size(-(sizeStimCm/2):0.5:(sizeStimCm/2))),-(sizeStimCm/2):0.5:(sizeStimCm/2));
                        worldPoints = table();
                        worldPoints.X = X(:);
                        worldPoints.Y = Y(:);
                        worldPoints.Z = zeros(size(worldPoints.X));



                    case 'CROSS'
                        % TODO: need to update the update funciton so the Z is not overwritten by a
                        % constant distance. Just shift things.
                        [X,Y] = meshgrid(-(sizeStimCm/2):0.5:(sizeStimCm/2),(0)*ones(size(-(sizeStimCm/2):0.5:(sizeStimCm/2))));
                        worldPoints = table();
                        worldPoints.X = X(:);
                        worldPoints.Y = Y(:);
                        worldPoints.Z = zeros(size(worldPoints.X));
                        % TODO: need to update the update funciton so the Z is not overwritten by a
                        % constant distance. Just shift things.
                        [X,Y] = meshgrid((0)*ones(size(-(sizeStimCm/2):0.5:(sizeStimCm/2))),-(sizeStimCm/2):0.5:(sizeStimCm/2));
                        worldPoints2 = table();
                        worldPoints2.X = X(:);
                        worldPoints2.Y = Y(:);
                        worldPoints2.Z = zeros(size(worldPoints2.X));

                        worldPoints = cat(1,worldPoints, worldPoints2);
                end
                app.Data.worldPoints = worldPoints;
                app.Data.stimulus = Values.Stimulus;
            else
                worldPoints = app.Data.worldPoints;
            end



            %% Update the points, eyes and screen acordint to the sliders
            app.Data.FixationSpot.X = Values.FixationX;
            app.Data.FixationSpot.Y = Values.FixationY;
            app.Data.FixationSpot.Z = Values.FixationDistance;
            app.Data.Screen.Slant = Values.ScreenSlant;

            leftEyeScreen = Geometry3D.MakeScreen(app.Data.Screen.SizeCm, app.Data.Screen.ResPix, app.Data.Screen.Distance, app.Data.Screen.Slant);
            rightEyeScreen = Geometry3D.MakeScreen(app.Data.Screen.SizeCm, app.Data.Screen.ResPix, app.Data.Screen.Distance, app.Data.Screen.Slant);

            eyes = Geometry3D.MakeEyes(Values.IPDMm/10, app.Data.FixationSpot, Values.TorsionVersion, Values.TorsionVergence);

            % Rotate the world points to apply the slant (rotates around
            % 0,0,0)
            worldPoints.X = worldPoints.X*Values.StimulusScale;
            worldPoints.Y = worldPoints.Y*Values.StimulusScale;
            worldPoints.Z = zeros(size(worldPoints.Z));

            % Rotate by slant and tilt
            R = Geometry3D.Quat2RotMat(Geometry3D.AxisAngle2Quat([cosd(Values.PlaneTilt) sind(Values.PlaneTilt) 0],deg2rad(Values.PlaneSlant)));
            worldPoints{:,:} = (R*worldPoints{:,:}')';

            % Displace by distance
            worldPoints.Z = worldPoints.Z + Values.StimulusDistance;

            % Add the fixation spot to the world points to conver to eye
            % and screen points.
            worldPoints{end+1,:} = [app.Data.FixationSpot.X app.Data.FixationSpot.Y app.Data.FixationSpot.Z];

            %% Get the points in eye and screen coordinates
            eyePoints = Geometry3D.Points3DToEyes(worldPoints, eyes);
            screenPoints = Geometry3D.PointsEyesToScreen(eyes, eyePoints, leftEyeScreen, rightEyeScreen);


            if ( ~isfield(app.Data, "f") || ~isvalid(app.Data.f))
                % If figure does not exist create it with all the plots and
                % handles to them
                [f, heyes, hfix, hscreen, hpoints, hspoints, hLRpoints, hdisparity, hscreen] = Geometry3D.demoDisparityInitPlots(worldPoints, eyePoints, screenPoints, eyes, leftEyeScreen, rightEyeScreen);
                app.Data.f = f;
                app.Data.heyes = heyes;
                app.Data.hfix = hfix;
                app.Data.hscreen = hscreen;
                app.Data.hpoints = hpoints;
                app.Data.hspoints = hspoints;
                app.Data.hLRpoints = hLRpoints;
                app.Data.hdisparity = hdisparity;
                app.Data.hscreen = hscreen;
            end

            % update 3D plot
            lxfar = eyes.L.X - 10000*eyes.L.RM(2,1);
            lyfar = eyes.L.Y - 10000*eyes.L.RM(3,1);
            lzfar = eyes.L.Z + 10000*eyes.L.RM(1,1);

            rxfar = eyes.R.X - 10000*eyes.R.RM(2,1);
            ryfar = eyes.R.Y - 10000*eyes.R.RM(3,1);
            rzfar = eyes.R.Z + 10000*eyes.R.RM(1,1);

            lrx = eyes.L.X - 10*eyes.L.RM(2,2);
            llx = eyes.L.X + 10*eyes.L.RM(2,2);
            lry = eyes.L.Y - 10*eyes.L.RM(3,2);
            lly = eyes.L.Y + 10*eyes.L.RM(3,2);
            lrz = eyes.L.Z + 10*eyes.L.RM(1,2);
            llz = eyes.L.Z - 10*eyes.L.RM(1,2);

            rrx = eyes.R.X - 10*eyes.R.RM(2,2);
            rlx = eyes.R.X + 10*eyes.R.RM(2,2);
            rry = eyes.R.Y - 10*eyes.R.RM(3,2);
            rly = eyes.R.Y + 10*eyes.R.RM(3,2);
            rrz = eyes.R.Z + 10*eyes.R.RM(1,2);
            rlz = eyes.R.Z - 10*eyes.R.RM(1,2);

            rbx = eyes.R.X - 10*eyes.R.RM(2,3);
            rtx = eyes.R.X + 10*eyes.R.RM(2,3);
            rby = eyes.R.Y - 10*eyes.R.RM(3,3);
            rty = eyes.R.Y + 10*eyes.R.RM(3,3);
            rbz = eyes.R.Z + 10*eyes.R.RM(1,3);
            rtz = eyes.R.Z - 10*eyes.R.RM(1,3);


            lbx = eyes.L.X - 10*eyes.L.RM(2,3);
            ltx = eyes.L.X + 10*eyes.L.RM(2,3);
            lby = eyes.L.Y - 10*eyes.L.RM(3,3);
            lty = eyes.L.Y + 10*eyes.L.RM(3,3);
            lbz = eyes.L.Z + 10*eyes.L.RM(1,3);
            ltz = eyes.L.Z - 10*eyes.L.RM(1,3);

            set(app.Data.hfix, 'xdata', Values.FixationX);
            set(app.Data.hfix, 'ydata', Values.FixationDistance);
            set(app.Data.hfix, 'zdata', Values.FixationY);
            set(app.Data.heyes.l, ...
                'xdata', [eyes.L.X lxfar ...
                nan lbx ltx nan lrx llx], ...
                'ydata', [eyes.L.Z lzfar  ...
                nan lbz ltz nan lrz llz], ...
                'zdata', [eyes.L.Y lyfar ...
                nan lby lty nan lry lly]);
            set(app.Data.heyes.r, ...
                'xdata', [ eyes.R.X rxfar ...
                nan rbx rtx nan rrx rlx ], ...
                'ydata', [ eyes.R.Z rzfar ...
                nan rbz rtz nan rrz rlz ], ...
                'zdata', [ eyes.R.Y ryfar ...
                nan rby rty nan rry rly]);
            set(app.Data.heyes.c(1), 'xdata', [eyes.L.X]);
            set(app.Data.heyes.c(2), 'xdata', [eyes.R.X]);

            set(app.Data.hpoints, 'xdata', worldPoints.X, 'ydata',worldPoints.Z, 'zdata', worldPoints.Y);
            set(app.Data.hspoints, 'xdata', worldPoints.X, 'ydata',worldPoints.Z, 'zdata', app.Data.hspoints.Parent.XLim(1)*ones(size(worldPoints.Y)));

            ax3d = app.Data.hfix.Parent;
            switch(app.Values.View3D)
                case 'Oblique'
                    view(ax3d, -60,10)
                case 'TOP'
                    view(ax3d, -90,90)
                case 'SIDE'
                    view(ax3d, -90,0)
            end

            % update retina plot
            ra = atan2( tand(eyePoints.RV) , tand(eyePoints.RH)./(abs(cosd(eyePoints.RV))));
            re = atand(sqrt(tand(eyePoints.RV).^2 + (tand(eyePoints.RH)./(abs(cosd(eyePoints.RV)))).^2));

            set(app.Data.hLRpoints.R, 'ThetaData', ra, 'RData', re);
            ra = atan2( tand(eyePoints.LV) , tand(eyePoints.LH)./abs(cosd(eyePoints.LV)));
            re = atand(sqrt(tand(eyePoints.LV).^2 + (tand(eyePoints.LH)./(abs(cosd(eyePoints.LV)))).^2));
            set(app.Data.hLRpoints.L, 'ThetaData', ra, 'RData', re);

            % update disparity plot
            set(app.Data.hdisparity, ...
                'xdata',eyePoints.RH, 'ydata', eyePoints.RV, ...
                'udata', eyePoints.LH-eyePoints.RH, 'vdata', eyePoints.LV-eyePoints.RV);

            % update screen plots
            set(app.Data.hscreen.LPoints, 'xdata', screenPoints.LX(1:end-1), 'ydata',screenPoints.LY(1:end-1));
            set(app.Data.hscreen.RPoints, 'xdata', screenPoints.RX(1:end-1), 'ydata',screenPoints.RY(1:end-1));
            set(app.Data.hscreen.LFP, 'xdata', screenPoints.LX(end), 'ydata', screenPoints.LY(end));
            set(app.Data.hscreen.RFP, 'xdata', screenPoints.RX(end), 'ydata', screenPoints.RY(end));

        end

    end


    methods(Static, Access = private) % DEMO MOTION FLOW

        function [f, heyes, hfix, hscreen, hpoints, hspoints, hLRpoints, hdisparity, hLscreen] = demoMotionFlowInitPlots(points, eyePoints, screenPoints, eyes, leftEyeScreen, rightScreen)

            t = points(1:end-1,:);
            fp = points(end,:);

            % View the dots, eyes, and fixation dot
            scr_siz = get(0,'ScreenSize');
            margin = floor(0.1*(scr_siz(4)));
            f = figure('color','w','position',floor([...
                margin...
                margin...
                scr_siz(3)*2.8/4 ...
                scr_siz(4)*2/4 ...
                ]));
            axes('OuterPosition',[ 0    0    0.5    1], 'nextplot','add')

            hpoints = line(0,0, 0,'linestyle','none','marker','o','Color','k');
            hspoints = line(0,0, 0,'linestyle','none','marker','o','Color',[0.8 0.8 0.8]);
            hfix = line(0,0, 0,'linestyle','none','marker','o','Color','r','LineWidth',2, 'markersize',20);

            heyes.c(1) = plot3([eyes.L.X], [eyes.L.Y ], 0,'o','Color','b', 'markersize',10); % left eye fixation spot and right eye
            heyes.c(2) = plot3([eyes.R.X], [eyes.R.Y], 0,'o','Color','r', 'markersize',10); % left eye fixation spot and right eye
            heyes.l = line(0,0,0,'linestyle','-','Color','b'); % left eye fixation spot and right eye
            heyes.r = line(0,0,0,'linestyle','-','Color','r'); % left eye fixation spot and right eye

            grid
            % xlim(1.2*[min([t.X;t.Y]) max([t.X;t.Y])])
            % zlim(1.2*[min([t.X;t.Y]) max([t.X;t.Y])])
            xlim([-100 100])
            zlim([-100 100])
            % ylim([-5 min(100,max(fp.Z, max(t.Z)))*2])
            ylim([-5 200])
            xlabel('X (cm)')
            zlabel('Y (cm)')
            ylabel('Z (cm)')
            view(-60,10)
            title('3D world');
            hscreen = [];
            % hscreen(1) = plot3([-1 -1 1 1 -1]*leftEyeScreen.widthCm/2+eyes.L.X, [1 1 1 1 1]*leftEyeScreen.ScreenDistance, [-1 1 1 -1 -1]*leftEyeScreen.heightCm/2);
            % hscreen(2) = plot3([-1 -1 1 1 -1]*rightScreen.widthCm/2+eyes.R.X, [1 1 1 1 1]*rightScreen.ScreenDistance, [-1 1 1 -1 -1]*rightScreen.heightCm/2);
            hLRpoints = [];
            t = eyePoints;

            % polaraxes('OuterPosition',[0.44    0.3    0.28    0.7], 'nextplot','add')
            % hLRpoints.L = polarplot(zeros(size(t.LV)), zeros(size(t.LV)),'o');
            % hLRpoints.R = polarplot(zeros(size(t.RV)), zeros(size(t.RV)),'o');
            % hLRpoints.FP = polarplot(0,0,'ro','linewidth',3);
            % grid
            % % xlim([min([t.RH t.LH],[],"all") max([t.RH t.LH],[],"all")])
            % % ylim([min([t.RV t.LV],[],"all") max([t.RV t.LV],[],"all")])
            % % xlabel('X (deg)')
            % % ylabel('Y (deg)')
            % set(gca,"RLim",[0 90])
            % set(gca,'ThetaTickLabel',[])
            % grid
            % % legend({'Right Eye','Left Eye'})
            % title('Visual field Left eye (blue) Right eye (red)');

            %-----------------------------
            % disparity plot
            %-----------------------------


            axes('OuterPosition',[ 0.65    0.1    0.35    0.9], 'nextplot','add')
            set(gca, 'PlotBoxAspectRatio',[1 1 1])
            hdisparity = quiver(t.RH, t.RV, t.LH-t.RH, t.LV-t.RV, 'AutoScale', "off");
            grid
            title('Disparity (Helmholtz, L->R)')
            set(gca,'xlim',[-40 40],'ylim',[-40 40])

        end

        function demoMotionFlowUpdate(app)

            Values = app.Values;

            if (~isfield(app.Data,"stimulus"))
                app.Data.stimulus = "NONE";
            end

            %%
            sizeStimCm = 40;

            %             ["Ground plane" "Cloud"])

            if (app.Data.stimulus ~= string(Values.Stimulus) )
                % Make some dots in world coordinates (XYZ)
                %
                % Make a table to start putting everything together
                switch (Values.Stimulus)
                    case 'Cloud'
                        numDots = 200;
                        worldPoints = table();
                        worldPoints.X = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
                        worldPoints.Y = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
                        worldPoints.Z = zeros(numDots, 1);
                    case 'Ground plane'
                        % TODO: need to update the update funciton so the Z is not overwritten by a
                        % constant distance. Just shift things.
                        [X,Y] = meshgrid([-(sizeStimCm/2):2:(sizeStimCm/2)],[-(sizeStimCm/2):2:(sizeStimCm/2)]);
                        worldPoints = table();
                        worldPoints.X = X(:);
                        worldPoints.Y = Y(:);
                        worldPoints.Z = zeros(size(worldPoints.X));
                end
                app.Data.worldPoints = worldPoints;
                app.Data.stimulus = Values.Stimulus;
            else
                worldPoints = app.Data.worldPoints;
            end



            %% Update the points, eyes and screen acordint to the sliders
            app.Data.FixationSpot.X = Values.FixationX;
            app.Data.FixationSpot.Y = Values.FixationY;
            app.Data.FixationSpot.Z = Values.FixationDistance;
            app.Data.Screen.Slant = Values.ScreenSlant;

            leftEyeScreen = Geometry3D.MakeScreen(app.Data.Screen.SizeCm, app.Data.Screen.ResPix, app.Data.Screen.Distance, app.Data.Screen.Slant);
            rightEyeScreen = Geometry3D.MakeScreen(app.Data.Screen.SizeCm, app.Data.Screen.ResPix, app.Data.Screen.Distance, app.Data.Screen.Slant);

            eyes = Geometry3D.MakeEyes(Values.IPDMm/10, app.Data.FixationSpot, Values.TorsionVersion, Values.TorsionVergence);

            % Rotate the world points to apply the slant (rotates around
            % 0,0,0)
            worldPoints.X = worldPoints.X*Values.StimulusScale;
            worldPoints.Y = worldPoints.Y*Values.StimulusScale;
            worldPoints.Z = zeros(size(worldPoints.Z));

            % Rotate by slant and tilt
            R = Geometry3D.Quat2RotMat(Geometry3D.AxisAngle2Quat([cosd(Values.PlaneTilt) sind(Values.PlaneTilt) 0],deg2rad(Values.PlaneSlant)));
            worldPoints{:,:} = (R*worldPoints{:,:}')';

            % Displace by distance
            worldPoints.Z = worldPoints.Z + Values.StimulusDistance;

            % Add the fixation spot to the world points to conver to eye
            % and screen points.
            worldPoints{end+1,:} = [app.Data.FixationSpot.X app.Data.FixationSpot.Y app.Data.FixationSpot.Z];

            %% Get the points in eye and screen coordinates
            eyePoints = Geometry3D.Points3DToEyes(worldPoints, eyes);
            screenPoints = Geometry3D.PointsEyesToScreen(eyes, eyePoints, leftEyeScreen, rightEyeScreen);


            if ( ~isfield(app.Data, "f") || ~isvalid(app.Data.f))
                % If figure does not exist create it with all the plots and
                % handles to them
                [f, heyes, hfix, hscreen, hpoints, hspoints, hLRpoints, hdisparity, hscreen] = Geometry3D.demoDisparityInitPlots(worldPoints, eyePoints, screenPoints, eyes, leftEyeScreen, rightEyeScreen);
                app.Data.f = f;
                app.Data.heyes = heyes;
                app.Data.hfix = hfix;
                app.Data.hscreen = hscreen;
                app.Data.hpoints = hpoints;
                app.Data.hspoints = hspoints;
                app.Data.hLRpoints = hLRpoints;
                app.Data.hdisparity = hdisparity;
                app.Data.hscreen = hscreen;
            end

            % update 3D plot
            lxfar = eyes.L.X - 10000*eyes.L.RM(2,1);
            lyfar = eyes.L.Y - 10000*eyes.L.RM(3,1);
            lzfar = eyes.L.Z + 10000*eyes.L.RM(1,1);

            rxfar = eyes.R.X - 10000*eyes.R.RM(2,1);
            ryfar = eyes.R.Y - 10000*eyes.R.RM(3,1);
            rzfar = eyes.R.Z + 10000*eyes.R.RM(1,1);

            lrx = eyes.L.X - 10*eyes.L.RM(2,2);
            llx = eyes.L.X + 10*eyes.L.RM(2,2);
            lry = eyes.L.Y - 10*eyes.L.RM(3,2);
            lly = eyes.L.Y + 10*eyes.L.RM(3,2);
            lrz = eyes.L.Z + 10*eyes.L.RM(1,2);
            llz = eyes.L.Z - 10*eyes.L.RM(1,2);

            rrx = eyes.R.X - 10*eyes.R.RM(2,2);
            rlx = eyes.R.X + 10*eyes.R.RM(2,2);
            rry = eyes.R.Y - 10*eyes.R.RM(3,2);
            rly = eyes.R.Y + 10*eyes.R.RM(3,2);
            rrz = eyes.R.Z + 10*eyes.R.RM(1,2);
            rlz = eyes.R.Z - 10*eyes.R.RM(1,2);

            rbx = eyes.R.X - 10*eyes.R.RM(2,3);
            rtx = eyes.R.X + 10*eyes.R.RM(2,3);
            rby = eyes.R.Y - 10*eyes.R.RM(3,3);
            rty = eyes.R.Y + 10*eyes.R.RM(3,3);
            rbz = eyes.R.Z + 10*eyes.R.RM(1,3);
            rtz = eyes.R.Z - 10*eyes.R.RM(1,3);


            lbx = eyes.L.X - 10*eyes.L.RM(2,3);
            ltx = eyes.L.X + 10*eyes.L.RM(2,3);
            lby = eyes.L.Y - 10*eyes.L.RM(3,3);
            lty = eyes.L.Y + 10*eyes.L.RM(3,3);
            lbz = eyes.L.Z + 10*eyes.L.RM(1,3);
            ltz = eyes.L.Z - 10*eyes.L.RM(1,3);

            set(app.Data.hfix, 'xdata', Values.FixationX);
            set(app.Data.hfix, 'ydata', Values.FixationDistance);
            set(app.Data.hfix, 'zdata', Values.FixationY);
            set(app.Data.heyes.l, ...
                'xdata', [eyes.L.X lxfar ...
                nan lbx ltx nan lrx llx], ...
                'ydata', [eyes.L.Z lzfar  ...
                nan lbz ltz nan lrz llz], ...
                'zdata', [eyes.L.Y lyfar ...
                nan lby lty nan lry lly]);
            set(app.Data.heyes.r, ...
                'xdata', [ eyes.R.X rxfar ...
                nan rbx rtx nan rrx rlx ], ...
                'ydata', [ eyes.R.Z rzfar ...
                nan rbz rtz nan rrz rlz ], ...
                'zdata', [ eyes.R.Y ryfar ...
                nan rby rty nan rry rly]);
            set(app.Data.heyes.c(1), 'xdata', [eyes.L.X]);
            set(app.Data.heyes.c(2), 'xdata', [eyes.R.X]);

            set(app.Data.hpoints, 'xdata', worldPoints.X, 'ydata',worldPoints.Z, 'zdata', worldPoints.Y);
            set(app.Data.hspoints, 'xdata', worldPoints.X, 'ydata',worldPoints.Z, 'zdata', app.Data.hspoints.Parent.XLim(1)*ones(size(worldPoints.Y)));

            ax3d = app.Data.hfix.Parent;
            switch(app.Values.View3D)
                case 'Oblique'
                    view(ax3d, -60,10)
                case 'TOP'
                    view(ax3d, -90,90)
                case 'SIDE'
                    view(ax3d, -90,0)
            end

            % update retina plot
            ra = atan2( tand(eyePoints.RV) , tand(eyePoints.RH)./(abs(cosd(eyePoints.RV))));
            re = atand(sqrt(tand(eyePoints.RV).^2 + (tand(eyePoints.RH)./(abs(cosd(eyePoints.RV)))).^2));

            set(app.Data.hLRpoints.R, 'ThetaData', ra, 'RData', re);
            ra = atan2( tand(eyePoints.LV) , tand(eyePoints.LH)./abs(cosd(eyePoints.LV)));
            re = atand(sqrt(tand(eyePoints.LV).^2 + (tand(eyePoints.LH)./(abs(cosd(eyePoints.LV)))).^2));
            set(app.Data.hLRpoints.L, 'ThetaData', ra, 'RData', re);

            % update disparity plot
            set(app.Data.hdisparity, ...
                'xdata',eyePoints.RH, 'ydata', eyePoints.RV, ...
                'udata', eyePoints.LH-eyePoints.RH, 'vdata', eyePoints.LV-eyePoints.RV);

            % update screen plots
            set(app.Data.hscreen.LPoints, 'xdata', screenPoints.LX(1:end-1), 'ydata',screenPoints.LY(1:end-1));
            set(app.Data.hscreen.RPoints, 'xdata', screenPoints.RX(1:end-1), 'ydata',screenPoints.RY(1:end-1));
            set(app.Data.hscreen.LFP, 'xdata', screenPoints.LX(end), 'ydata', screenPoints.LY(end));
            set(app.Data.hscreen.RFP, 'xdata', screenPoints.RX(end), 'ydata', screenPoints.RY(end));

        end

    end

    methods(Static, Access = private) % DEMO QUATERNIONS

        function [f, h] = demoQuaternionInitPlots(app)
            scr_siz = get(0,'ScreenSize');
            margin = floor(0.1*(scr_siz(4)));
            f = figure('color','w','position',floor([...
                margin...
                margin...
                scr_siz(3)*2/4 ...
                scr_siz(4)*2/4 ...
                ]));

            % Define the radius of the sphere
            R = 1;
            % Define the resolution of the sphere
            res = 20;
            % Define the transparency of the sphere
            alpha = 0.8;

            % Create the sphere
            [Z,X,Y] = sphere(res);

            % Plot the sphere
            h.sphere = mesh(R*X,R*Y,R*Z,'FaceAlpha', alpha,'facecolor',[0.8 0.8 0.8],'edgecolor',0.2*[0.8 0.8 0.8]);
            axis equal;

            %
            % set(gcf, 'Renderer', 'OpenGL');
            % shading interp, material shiny, lighting phong, lightangle(0, 55);

            view(45,30)
            set(gca,'xlim',1.5*[-1 1],'ylim',1.5*[-1 1],'zlim',1.5*[-1 1])

            h.ax = gca;

            h.frame(1) = line(0,0,0,'linewidth',2,'color','k');
            h.frame(2) = line(0,0,0,'linewidth',2,'color','r');
            h.frame(3) = line(0,0,0,'linewidth',2,'color','k');
            h.rotax = line(0,0,0,'linestyle','-','linewidth',2,'color','b');

        end

        function demoQuaternionUpdate(app)

            Values = app.Values;
            q = [Values.q0 Values.q1 Values.q2 Values.q3];
            hvt = [-Values.H -Values.V Values.T];

            if ( ~isfield(app.Data, "f") || ~isvalid(app.Data.f))
                % If figure does not exist create it with all the plots and
                % handles to them
                [f, h] = Geometry3D.demoQuaternionInitPlots();
                app.Data.f = f;
                app.Data.h = h;

                app.Data.LastQ = q;
                app.Data.Lasthvt = hvt;
                app.Data.LastCS = app.Values.x2DHVTCoordinateSystem;
            end



            steps = [0:10:360];
            points = -90:10:90;

            [az, el] = meshgrid(steps,points);


            q1 = app.Data.LastQ;

            % Update only the 3 values not changed by the user
            % so as to keep a valid quaternion. If the other three values
            % are zero transform to aligned with 1 0 0 because otherwise
            % there is a singularity
            c = (q-q1) == 0;
            %sum(c)
            if ( sum(c)==3)
                nc = q(c);
                if ( sum(nc) == 0)
                    nc = [1 0 0];
                end

                q(c) = nc*sqrt(1-q(~c)^2)/norm(nc);

                switch(app.Values.x2DHVTCoordinateSystem)
                    case 'Fick'
                        hvt = rad2deg(Geometry3D.RotMat2Fick(Geometry3D.Quat2RotMat(q)));
                    case 'Helmholtz'
                        hvt = rad2deg(Geometry3D.RotMat2Helm(Geometry3D.Quat2RotMat(q)));
                    case 'Listings(AET)'
                        hvt = rad2deg(Geometry3D.RotMat2Listings(Geometry3D.Quat2RotMat(q)));
                end


            elseif (~isequal(hvt, app.Data.Lasthvt)  )
                % if the quaternion was not update check if the hvt was

                switch(app.Values.x2DHVTCoordinateSystem)
                    case 'Fick'
                        q = Geometry3D.RotMat2Quat(Geometry3D.Fick2RotMat(deg2rad(hvt)));
                    case 'Helmholtz'
                        q = Geometry3D.RotMat2Quat(Geometry3D.Helm2RotMat(deg2rad(hvt)));
                    case 'Listings(AET)'
                        q = Geometry3D.RotMat2Quat(Geometry3D.Listings2RotMat(deg2rad(hvt)));
                end
            elseif (string(app.Data.LastCS) ~= string(app.Values.x2DHVTCoordinateSystem))

                switch(app.Values.x2DHVTCoordinateSystem)
                    case 'Fick'
                        hvt = rad2deg(Geometry3D.RotMat2Fick(Geometry3D.Quat2RotMat(q)));
                    case 'Helmholtz'
                        hvt = rad2deg(Geometry3D.RotMat2Helm(Geometry3D.Quat2RotMat(q)));
                    case 'Listings(AET)'
                        hvt = rad2deg(Geometry3D.RotMat2Listings(Geometry3D.Quat2RotMat(q)));
                end
            end

            xlabel(app.Data.h.ax,'x')
            ylabel(app.Data.h.ax,'y')
            zlabel(app.Data.h.ax,'z')

            switch(app.Values.x2DHVTCoordinateSystem)
                case 'Fick'
                    z = sind(az);
                    y = sind(el).*cosd(az);
                    x = cosd(az).*cosd(el);
                case 'Helmholtz'
                    y = sind(az);
                    x = sind(el).*cosd(az);
                    z = cosd(az).*cosd(el);
                case 'Listings(AET)'
                    x = sind(az);
                    z = sind(el).*cosd(az);
                    y = cosd(az).*cosd(el);

                case 'Harms'
                    x = sqrt(1./(1+tand(az).^2+tand(el).^2));
                    x(az > 90 | az < -90) = -x(az > 90 | az < -90);
                    x(az == 90 | az == -90 | el == 90 | el == -90) = 0;

                    z = tand(el).*x;
                    y = tand(az).*x;

                    z(el == 90 | el == -90) = 1;
                    z(el == 90 | el == -90) = 1;
                    y(el == 90 | el == -90) = 0;

                    z(az == 90 | az == -90) = 0;
                    y(az == 90 | az == -90) = 1;
                case 'Hess'
                    z = sind(el);
                    y = sind(az);
                    x = sqrt(1-z.^2-y.^2);
                    x(az >= 90 | az <= -90) = -x(az >= 90 | az <= -90);
                    x(az == 90 | az == -90) = 0;

                    outsidepoints = (z.^2+y.^2)>=1;
                    z(outsidepoints) = nan;
                    y(outsidepoints) = nan;
                    x(outsidepoints) = nan;
                    % x = real(x);
            end

            set(app.Data.h.sphere, 'xdata',x, 'ydata',y, 'zdata',z)
            %             set(app.Data.h)
            %

            app.Data.Lasthvt = hvt;
            app.Data.LastQ = q;
            app.Data.LastCS = app.Values.x2DHVTCoordinateSystem;

            app.Values.q0  = q(1);
            app.Values.q1  = q(2);
            app.Values.q2  = q(3);
            app.Values.q3  = q(4);

            app.Values.H  = -hvt(1);
            app.Values.V  = -hvt(2);
            app.Values.T  = hvt(3);
            % app.Values.x  = axis(1);
            % app.Values.y  = axis(2);
            % app.Values.z  = axis(3);
            % app.Values.a  = angle;


            R = Geometry3D.Quat2RotMat(q);
            c = R(:,1);
            axis = R(:,3);
            %
            set(app.Data.h.frame(1), 'xdata',1.2*[0 R(1,1)], 'ydata',1.2*[0 R(2,1)], 'zdata',1.2*[0 R(3,1)])
            set(app.Data.h.frame(2), 'xdata',c(1)+ 0.5*[0 R(1,2)], 'ydata',c(2)+ 0.5*[0 R(2,2)], 'zdata',c(3)+ 0.5*[0 R(3,2)])
            set(app.Data.h.frame(3), 'xdata',c(1)+ 0.2*[0 R(1,3)], 'ydata',c(2)+ 0.2*[0 R(2,3)], 'zdata',c(3)+ 0.2*[0 R(3,3)])

            set(app.Data.h.rotax, 'xdata', 1.2*axis(1)*[0 1], 'ydata', 1.2*axis(2)*[0 1], 'zdata', 1.2*axis(3)*[0 1])
            %
            % h.frame(1) = line(0,0,0,'linewidth',2);
            % h.frame(2) = line(0,0,0,'linewidth',2);
            % h.frame(3) = line(0,0,0,'linewidth',2);


        end
    end

    methods(Static)
        % rotations around head fixed axis
        % S = [sx, sy, sz] is a right handed space fixed coordinate system
        %    sx line of sight (pointing forward)
        %    sy interarual axis (pointing left)
        %    sz earth vertical (pointing up)
        %
        % B = [bx, by, bz] is a right handed eye fixed coordinate system
        %
        %

        % HORIZONTAL ROTATION right handed
        function M = RotZ(theta)
            M = [   cos(theta)  -sin(theta)     0;
                sin(theta)  cos(theta)      0;
                0           0               1];
        end

        % VERTICAL ROTATION right handed
        function M = RotY(phi)
            M = [   cos(phi)  0               sin(phi);
                0           1               0;
                -sin(phi) 0               cos(phi)];
        end

        % TORSIONAL ROTATION right handed
        function M = RotX(psi)
            M = [   1           0               0;
                0           cos(psi)      -sin(psi);
                0           sin(psi)      cos(psi)];
        end
    end

    methods(Static)

        % Fick to rotation matrix
        function M = Fick2RotMat(HVT)

            H = HVT(1);
            V = HVT(2);
            T = HVT(3);

            M = Geometry3D.RotZ(H)*Geometry3D.RotY(V)*Geometry3D.RotX(T);
        end

        % Rotation matrix to Fick
        function HVT = RotMat2Fick(M)
            % TODO: Double check this.
            r31 = M(3,1);
            r21 = M(2,1);
            r32 = M(3,2);

            HVT(2) = -asin(r31);
            HVT(1) = asin(r21/cos(HVT(2)));
            HVT(3) = asin(r32/cos(HVT(2)));
        end

        function M = Helm2RotMat(HVT)

            H = HVT(1);
            V = HVT(2);
            T = HVT(3);

            M = Geometry3D.RotY(V)*Geometry3D.RotZ(H)*Geometry3D.RotX(T);
        end

        function HVT = RotMat2Helm(M)
            r21 = M(2,1);
            r31 = M(3,1);
            r23 = M(2,3);

            HVT(1) = asin(r21);
            HVT(2) = -asin(r31/cos(HVT(1)));
            HVT(3) = -asin(r23/cos(HVT(1)));
        end

        function M = List2Mat(ADT)
            A = ADT(1);
            D = ADT(2);
            T = ADT(3);
            M = Geometry3D.RotX(A)*Geometry3D.RotZ(D)*Geometry3D.RotX(T-A);
        end

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

        function M = Quat2RotMat(q)
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

            M = E*G';

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

        %
        %         function M = Quat2RotMat(q)
        %         end
        %
        %         function HVT = RotMat2Fick(M)
        %         end
        %
        %         function HVT = RotMat2Helm(M)
        %         end
        %
        %         function q = Mat2


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
            x = a(a(:,1)>0,1);
            y = a(a(:,1)>0,2);
            z = a(a(:,1)>0,3);

            [dudx, dudy, dudz, dvdx, dvdy, dvdz] = Geometry3D.TangentSphereLinearJacobian(x,y,z);
            [duwxdt, duwydt, duwzdt, dvwxdt, dvwydt, dvwzdt] = Geometry3D.TangentSphereRotationalJacobian(x,y,z);

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
            hl(1) = quiver3(x,y,z, dudx, dudy, dudz,'linewidth',2,'DisplayName', '\partial xyz / \partial u');
            hl(2) = quiver3(x,y,z, dvdx, dvdy, dvdz,'linewidth',2,'DisplayName', '\partial xyz / \partial v');
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
            [x,y,z] = Geometry3D.FickToSphere(deg2rad([0 50]),deg2rad([0 40]));


            [dudx, dudy, dudz, dvdx, dvdy, dvdz] = Geometry3D.TangentSphereLinearJacobian(x,y,z);

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
            % quiver3(x,y,z, x, y, z,'linewidth',1)

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
                type = 'NUMERICAL';
            end

            step = 10;
            range = 80;
            [az, el] = meshgrid(-range:step:range,-range:step:range); % azimuths and elevations to include
            coordinateSystems = {'sphere', 'Fick', 'Helmholtz', 'Harms','Hess'};

            f = figure('color','white');
            tiledlayout(5,7,"TileSpacing","tight","Padding","tight");


            f.Position = [f.Position(1) f.Position(2) 4*f.Position(3) f.Position(4)];

            for i=1:length(coordinateSystems)
                sys = coordinateSystems{i};

                nexttile

                set(gca, 'nextplot','add')
                % subplot(4,length(coordinateSystems),i,'nextplot','add');

                % calculate spherical coordinates depending on the coordinate system
                switch(sys)
                    case 'sphere'
                        v = Geometry3D.SampleVisualDirections(100);
                        x = v(:,1); y=v(:,2); z=v(:,3);
                    case 'Fick'
                        [x,y,z] = Geometry3D.FickToSphere(deg2rad(az),deg2rad(el));
                    case 'Helmholtz'
                        [x,y,z] = Geometry3D.HelmholtzToSphere(deg2rad(az),deg2rad(el));
                    case 'Harms'
                        [x,y,z] = Geometry3D.HarmsToSphere(deg2rad(az),deg2rad(el));
                    case 'Hess'
                        [x,y,z] = Geometry3D.HessToSphere(deg2rad(az),deg2rad(el));
                end

                if ( i>1)
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
                end


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
                dw = deg2rad(0.1); % differential angular rotation in deg

                az2 = az;
                el2 = el;

                colorsflow = {'r' 'b' 'k'};
                dims = {'x' 'y' 'z'};

                dxyz = {[dv,0,0],[0,dv,0],[0,0,dv]}; % assumes distance to the target equal to the radius. It should be scaled by the ratio between the depth and R
                dRM = {Geometry3D.RotX(dw), Geometry3D.RotY(dw), Geometry3D.RotZ(dw)};

                for iFlow = 1:2
                    switch(iFlow)
                        case 1
                            WHICHFLOW = 'linear';
                        case 2
                            WHICHFLOW = 'rotational';
                    end
                    for vi=1:3

                        switch(type)
                            case 'NUMERICAL'
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
                                        [az2, el2 ] = Geometry3D.SphereToFick(x2,y2,z2);
                                    case 'Helmholtz'
                                        [az2, el2 ] = Geometry3D.SphereToHelmholtz(x2,y2,z2);
                                    case 'Harms'
                                        [az2, el2 ] = Geometry3D.SphereToHarms(x2,y2,z2);
                                    case 'Hess'
                                        [az2, el2 ] = Geometry3D.SphereToHess(x2,y2,z2);
                                end

                                daz = rad2deg(az2)-az;
                                del = rad2deg(el2)-el;

                            case 'ANALYTICAL'

                                % calculate spherical coordinates depending on the coordinate system
                                switch(sys)
                                    case 'sphere'
                                        [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.TangentSphereLinearJacobian(x,y,z);
                                        [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = Geometry3D.TangentSphereRotationalJacobian(x,y,z);
                                    case 'Fick'
                                        [x,y,z] = Geometry3D.FickToSphere(deg2rad(az),deg2rad(el));
                                        [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.FickLinearJacobian(deg2rad(az), deg2rad(el));
                                        [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = Geometry3D.FickRotationalJacobian(deg2rad(az), deg2rad(el));
                                    case 'Helmholtz'
                                        [x,y,z] = Geometry3D.HelmholtzToSphere(deg2rad(az),deg2rad(el));
                                        [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HelmholtzLinearJacobian(deg2rad(az), deg2rad(el));
                                        [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = Geometry3D.HelmholtzRotationalJacobian(deg2rad(az), deg2rad(el));

                                    case 'Harms'
                                        [x,y,z] = Geometry3D.HarmsToSphere(deg2rad(az),deg2rad(el));
                                        [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HarmsLinearJacobian(deg2rad(az), deg2rad(el));
                                        [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = Geometry3D.HarmsRotationalJacobian(deg2rad(az), deg2rad(el));
                                    case 'Hess'
                                        [x,y,z] = Geometry3D.HessToSphere(deg2rad(az),deg2rad(el));
                                        [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HessLinearJacobian(deg2rad(az), deg2rad(el));
                                        [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = Geometry3D.HessRotationalJacobian(deg2rad(az), deg2rad(el));
                                end

                                dazlin = {dazdx, dazdy, dazdz};
                                dellin = {deldx, deldy, deldz};
                                dazrot = {dazwxdt, dazwydt, dazwzdt};
                                delrot = {delwxdt, delwydt, delwzdt};

                                if ( strcmp( WHICHFLOW , 'linear' ) )
                                    daz = dazlin{vi};
                                    del = dellin{vi};
                                else
                                    daz = dazrot{vi};
                                    del = delrot{vi};
                                end

                                x2 = x;
                                y2 = y;
                                z2 = z;
                        end


                        nexttile
                        set(gca,'nextplot','add')
                        if ( i>1 )
                            quiver(az, el, rad2deg(daz), rad2deg(del), colorsflow{vi}, 'linewidth',1.5)

                            plot(0,0,'ro')

                            % if ( i > 1)
                            %     xlabel('Azimuth (deg)')
                            %     ylabel('Elevation (deg)')
                            % else
                            %     xlabel('Angle (deg)')
                            %     ylabel('Eccentricity (deg)')
                            % end
                            grid on

                            axis equal;
                            % set(gca,'xlim',[-82 82],'ylim',[-82 82])
                            set(gca,'xlim',[-92 92],'ylim',[-92 92])
                            set(gca,'xtick',[-90:30:90],'ytick',[-90:30:90])
                        else
                            view(140,15)
                            axis equal;

                            R = 1; % radius of the eye
                            step = 10;
                            range = 90;
                            [azs, els] = meshgrid(-range:step:range,-range:step:range); % azimuths and elevations to include
                            coordinateSystems = {'sphere', 'Fick', 'Helmholtz', 'Harms','Hess'};

                            [xs,ys,zs] = Geometry3D.FickToSphere(deg2rad(azs),deg2rad(els));
                            mesh(xs,ys,zs,'FaceAlpha', 0.5,'facecolor',0.8*[1 1 1],'EdgeColor','none');
                            % mesh(x,y,z,'FaceAlpha', 0.9,'facecolor',[1 1 1]);
                            % xlabel(gca,'x')
                            % ylabel(gca,'y')
                            % zlabel(gca,'z')
                            % line of sight
                            % line([0 R*1.5 ],[0 0 ],[0 0 ],'color','r','linewidth',2)
                            % title(sys)
                            quiver3(x,y,z,x2-x,y2-y,z2-z, 'linewidth',1.5)
                            
                            m(:,1,1) = [1 0 0 0.3 0 0];
                            m(:,1,2) = [0 1 0 0 0.3 0];
                            m(:,1,3) = [0 0 1 0 0 .3];
                            m(:,2,1) = [-2 0 0 2 0 0];
                            m(:,2,2) = [0 -2 0 0 2 0];
                            m(:,2,3) = [0 0 -2 0 0 2];
                            if ( iFlow==1)
                                quiver3(m(1,iFlow,vi),m(2,iFlow,vi),m(3,iFlow,vi),m(4,iFlow,vi),m(5,iFlow,vi),m(6,iFlow,vi),'color','r', 'linewidth',2)
                            else
                                 line([m(1,iFlow,vi) m(4,iFlow,vi)],[m(2,iFlow,vi) m(5,iFlow,vi)],[m(3,iFlow,vi) m(6,iFlow,vi)],'color','r','linewidth',2,'linestyle','-.')
                            end
                            set(gca,'xlim',[-1.3 1.3],'ylim',[-1.3 1.3],'zlim',[-1.3 1.3])

                            switch(WHICHFLOW)
                                case 'linear'
                                    title(['Flow due to ' WHICHFLOW  ' motion along ',dims{vi}])
                                case 'rotational'
                                    title(['Flow due to ' WHICHFLOW  ' motion around ',dims{vi}])
                            end
                        end
                    end
                end
            end
            delete(nexttile(1))
        end


        % Functions to converte between spherical coordinates and coordinate
        % sysstems for 2D rotations

        % Listings is a polar system with angle and eccentricity
        % Fick is a azimuth as latitudes (parallels) and elevation as longitudes (meridians)
        % Helmoltz is a azimuth as longitudes (meridians) and elevation as latitudes (parallels)
        % Harms is a azimuth as longitudes (meridians) and elevation as longitudes (meridians)
        % Hess is a azimuth as latitudes (parallels) and elevation as latitudes (parallels)

        function [x, y, z] = FickToSphere(az, el)
            x = cos( el ) .* cos( az );
            y = sin( az ) .* cos( el );   % longitude
            z = sin( el );                 % latitude
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

        function [x, y, z] = HarmsToSphere(az, el)
            azdeg = rad2deg(az);
            eldeg = rad2deg(el);
            x = 1 ./ sqrt(1 + tan(az).^2 + tan(el).^2 );
            x(azdeg > 90 | azdeg < -90) = -x(azdeg > 90 | azdeg < -90);
            x(azdeg == 90 | azdeg == -90 | eldeg == 90 | eldeg == -90) = 0;

            y = tan(az) .* x;
            z = tan(el) .* x;

            z(eldeg == 90 ) = 1;
            z( eldeg == -90) = 1;
            y(eldeg == 90 | eldeg == -90) = 0;

            z(azdeg == 90 | azdeg == -90) = 0;
            y(azdeg == 90 ) = 1;
            y(azdeg == -90) = -1;
        end

        function [az,el] = SphereToHarms(x,y,z)
            az = atan2(y,x);  % longitude
            el = atan2(z,x);  % longitude
        end


        function [dazdx, dazdy, dazdz, deldx, deldy, deldz] = FickLinearJacobian(az, el)
            % units should be deg/m

            [x, y, z] = Geometry3D.FickToSphere(az, el);
            % this are just the derivatives of the SphereToHarms function

            % D = sqrt(x.^2+y.^2+z.^2);
            % az = atan2d(y,x);
            % el = asind(z/D);

            dazdx =  -y./(y.^2+x.^2);
            dazdy = x./(y.^2+x.^2);
            dazdz = zeros(size(az));

            DD = sqrt(x.^2+y.^2);
            deldx = -z.*x ./ DD;
            deldy = -z.*y ./ DD;
            deldz =  DD;

        end

        function [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = FickRotationalJacobian(az, el)
            % unitless or deg per deg or rad per rad
            [x, y, z] = Geometry3D.FickToSphere(az, el);

            [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.FickLinearJacobian(az, el);

            % this uses the property that the linear velocity of a point on the
            % sphere rotating by an angular velocity is the cross product of the
            % coordinates of the point and the angular velocity.
            %
            % so then we can multiply the linear jacobian by the cross product
            % matrix equivalent of the point coordinates

            dazwxdt = 0*dazdx + z.*dazdy - y.*dazdz;
            dazwydt = -z.*dazdx - 0.*dazdy + x.*dazdz;
            dazwzdt = y.*dazdx - x.*dazdy + 0*dazdz;

            delwxdt = 0*deldx + z.*deldy - y.*deldz ;
            delwydt = -z.*deldx + 0*deldy + x.*deldz;
            delwzdt = y.*deldx - x.*deldy + 0.*deldz ;
        end

        function [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel]  = FickLinearInverseJacobian(az, el)

            dxdaz = -sin( az ) .* cos( el );
            dydaz = cos( az ) .* cos( el );
            dzdaz = 0;
            dxdel = -cos( az ) .* sin( el );
            dydel = -sin( az ) .* sin( el );
            dzdel = cos( el );
        end

        function [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HelmholtzLinearJacobian(az, el)
            % units should be deg/m

            [x, y, z] = Geometry3D.HelmholtzToSphere(az, el);
            % this are just the derivatives of the SphereToHarms function

            % D = sqrt(x.^2+y.^2+z.^2);
            % az = atan2d(y,x);
            % el = asind(z/D);

            DD = sqrt(x.^2+z.^2);
            dazdx = -y.*x ./ DD;
            dazdy =  DD;
            dazdz = -z.*y ./ DD;

            deldx =  -z./(z.^2+x.^2);
            deldy = zeros(size(az));
            deldz = x./(z.^2+x.^2);
        end

        function [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = HelmholtzRotationalJacobian(az, el)
            % unitless or deg per deg or rad per rad
            [x, y, z] = Geometry3D.HelmholtzToSphere(az, el);

            [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HelmholtzLinearJacobian(az, el);

            % this uses the property that the linear velocity of a point on the
            % sphere rotating by an angular velocity is the cross product of the
            % coordinates of the point and the angular velocity.
            %
            % so then we can multiply the linear jacobian by the cross product
            % matrix equivalent of the point coordinates

            dazwxdt = 0*dazdx + z.*dazdy - y.*dazdz;
            dazwydt = -z.*dazdx - 0.*dazdy + x.*dazdz;
            dazwzdt = y.*dazdx - x.*dazdy + 0*dazdz;

            delwxdt = 0*deldx + z.*deldy - y.*deldz ;
            delwydt = -z.*deldx + 0*deldy + x.*deldz;
            delwzdt = y.*deldx - x.*deldy + 0.*deldz ;
        end

        function [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel]  = HelmholtzLinearInverseJacobian(az, el)

            dxdaz = -sin(az) .* cos(el);
            dydaz = cos(az);
            dzdaz = -sin(el) .* sin(az);

            dxdel = -cos(az) .* sin(el);
            dydel = zeros(size(az));
            dzdel = cos(az) .* cos(el);

        end



        function [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HarmsLinearJacobian(az, el)
            [x, y, z] = Geometry3D.HarmsToSphere(az, el);

            % az = atan2d(y,x);
            % el = atan2d(z,x);
            %
            % x = 1/sqrt(1+tan(az)^2 +tan(el)^2)
            % y = tan(az)*x
            % z = tan(el)*x
            %
            %
            % this are just the derivatives of the SphereToHarms function
            dazdx = -y./(y.^2+x.^2);
            dazdy = x./(y.^2+x.^2);
            dazdz = zeros(size(az));

            deldx = -z./(z.^2+x.^2);
            deldy = zeros(size(el));
            deldz = x./(z.^2+x.^2);

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

        function [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = HarmsLinearInverseJacobian(az, el)
            % Define the trigonometric functions
            tan_az = tan(az);
            sec_az = sec(az);
            tan_el = tan(el);
            sec_el = sec(el);

            % Common denominator
            denominator = (1 + tan_az.^2 + tan_el.^2).^(3/2);

            % Calculate each element of the Jacobian matrix
            dxdaz = -tan_az .* sec_az.^2 ./ denominator;
            dxdel = -tan_el .* sec_el.^2 ./ denominator;

            dydaz = (sec_az.^2 - tan_az.^2 .* sec_az.^2) ./ denominator;
            dydel = -tan_az .* sec_az.^2 .* tan_el ./ denominator;

            dzdaz = -tan_az .* sec_az.^2 .* tan_el ./ denominator;
            dzdel = (sec_el.^2 - tan_el.^2 .* sec_el.^2) ./ denominator;

        end

        function [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = HarmsRotationalJacobian(az, el)
            [x, y, z] = Geometry3D.HarmsToSphere(az, el);

            [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HarmsLinearJacobian(az, el);

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

        function [dazdx, dazdy, dazdz, deldx, deldy, deldz] = HessLinearJacobian(az, el)
            [x, y, z] = Geometry3D.HessToSphere(az, el);
            % this are just the derivatives of the SphereToHarms function
            dazdx =  -x.*y./(sqrt(1-y.^2)) ;
            dazdy =  (x.^2 + z.^2 ) ./(sqrt(1-y.^2)) ;
            dazdz = -y.*z./(sqrt(1-y.^2)) ;

            deldx =  -x.*z./(sqrt(1-y.^2)) ;
            deldy = -y.*z./(sqrt(1-y.^2)) ;
            deldz =  (x.^2 + y.^2 ) ./(sqrt(1-y.^2)) ;
        end


        function [dazwxdt, dazwydt, dazwzdt, delwxdt, delwydt, delwzdt] = HessRotationalJacobian(az, el)
            [x, y, z] = Geometry3D.HessToSphere(az, el);

            [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HessLinearJacobian(az, el);

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



        function [x, y, z] = HessToSphere(az, el)
            azdeg = rad2deg(az);

            z = sin(el);
            y = sin(az);
            x = sqrt(1-z.^2-y.^2);
            % Need to flip negative azimuts
            x(azdeg >= 90 | azdeg <= -90) = -x(azdeg >= 90 | azdeg <= -90);
            % Need to force zero at 90 deg azimuth to avoid some complex
            % numbers that can appear numerically
            x(azdeg == 90 | azdeg == -90) = 0;

            % some of the combinations of azimuth and elevation actually
            % don't exist within the sphere.
            outsidepoints = (z.^2+y.^2)>=1;
            z(outsidepoints) = nan;
            y(outsidepoints) = nan;
            x(outsidepoints) = nan;
            % x = real(x);
        end

        function [az,el] = SphereToHess(x,y,z)
            D = sqrt( x.^2 + y.^2 + z.^2 );

            el = asin( z ./ D ); % latitude
            az = asin( y ./ D ); % latitude
        end

        function [x, y, z] = ListingsToSphere(angle, eccentricity)
            x = cos( eccentricity );
            y = sin( eccentricity ) .* cos( angle );
            z = sin( eccentricity ) .* sin( angle );
        end

        function [angle, ecc] = SphereToListings(x,y,z)
            ecc = acos(x);
            angle = atan2(z,y);
        end

        function [dudx, dudy, dudz, dvdx, dvdy, dvdz] = TangentSphereLinearJacobian(x,y,z)

            dudx = -y;
            dudy = 1 -  ( y.^2 ) ./ (1 + x);
            dudz = -( y .* z ) ./ (1 + x) ;

            dvdx = -z ;
            dvdy = -( y .* z ) ./ (1 + x) ;
            dvdz = 1 -  ( z.^2 ) ./ (1 + x) ;

        end

        function [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = TangentSphereInverseJacobian(x,y,z)
            [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.TangentSphereLinearJacobian(x,y,z);
        end

        function [duwxdt, duwydt, duwzdt, dvwxdt, dvwydt, dvwzdt] = TangentSphereRotationalJacobian(x,y,z)

            duwxdt = y;
            duwydt = 1 -  ( y.^2 ) ./ (1 + x);
            duwzdt = -( y .* z ) ./ (1 + x) ;

            dvwxdt = -z;
            dvwydt = -( y .* z ) ./ (1 + x);
            dvwzdt =  1 -  ( z.^2 ) ./ (1 + x);
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
end
% TODO:
% https://work.thaslwanter.at/thLib/html/rotmat.html
