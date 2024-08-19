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
            app.AddDropDown('Coordinate System',   1,  ["Fick", "Helmholtz", "Harms","Hess"])
            % app.AddSlider('Eye Radius',           0.02,  [0.01    1])
            %app.AddDropDown('Stimulus',      1,  ["Ground plane" "Point cloud"])
            app.AddSlider('Azimuth',           -40,  [-90 90])
            app.AddSlider('Elevation',           20,  [-90 90])
            % app.AddSlider('Torsion Version',      0,  [-20  20])
            % app.AddSlider('Torsion Vergence',     0,  [-20  20])
            %app.AddSlider('Ground plane slant',          0,  [-90  90])
            % app.AddSlider('Ground plane tilt',           0,  [0    90])
            app.AddSlider('Angular velocity X (deg/s)',   1*60,  [-100 100] )
            app.AddSlider('Angular velocity Y (deg/s)',   -0.5*60,  [-100 100])
            app.AddSlider('Angular velocity Z (deg/s)',   1.3*60,  [-100 100])
            app.AddSlider('Linear velocity X (m/s)',    .5,  [-5 5])
            app.AddSlider('Linear velocity Y (m/s)',    1,  [-5 5])
            app.AddSlider('Linear velocity Z (m/s)',    -1,  [-5 5])
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
        %%
        function demoCoordinateSystemsAndPlaneUpdate(app)


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
                app.Data.hs.textpoint = text(0*1.2,0*1.2,0*1.2, '(\theta,\psi)','fontsize',14, 'FontWeight','normal','HorizontalAlignment','right','VerticalAlignment','top');

                % draw tangent
                app.Data.hs.meshtangent = mesh(zeros(2,2), zeros(2,2),zeros(2,2), 'EdgeColor',[0.5 0.5 0.5],'FaceAlpha', 0.3,'facecolor',[0.8 0.8 0.8]);
                app.Data.hs.quiverdaz = quiver3(0,0,0, 0*2, 0*2 ,0*2,'color','k','linewidth',2);
                app.Data.hs.quiverdel = quiver3(0,0,0, 0*2, 0*2 ,0*2,'color','k','linewidth',2 );
                app.Data.hs.quivertvw = quiver3(0,0,0, 0*2, 0*2 ,0*2,'color',colors(3,:),'linewidth',2  ,'LineStyle',':');
                app.Data.hs.quivertvv = quiver3(0,0,0, 0*2, 0*2 ,0*2,'color',colors(4,:),'linewidth',2  ,'LineStyle',':');
                app.Data.hs.quivertv  = quiver3(0,0,0, 0*2, 0*2 ,0*2,'color','k','linewidth',2 ,'LineStyle','--');


                % draw flat
                app.Data.hs.ax2 = subplot(1,2,2,'nextplot','add');
                axis equal;

%                 app.Data.hs.meshflat = mesh(zeros(2,2), zeros(2,2), zeros(2,2) ,'FaceAlpha', 0.9,'facecolor',[1 1 1]);

                set(gca,'xtick',[],'ytick',[])
%                 view(90,0)
                %                 set(gca,'visible','off')
%                 title(altTitles{i})
                xlabel('Azimuth ')
                ylabel('Elevantion ')

                set(gca,'xlim',[-1.5 1.5],'ylim',[-1.5 1.5])

%                 app.Data.hs.meshpointflat = mesh(zeros(2,2), zeros(2,2),zeros(2,2), 'EdgeColor','none','FaceColor','k');
                app.Data.hs.textpointflat = text(0*1.2,0*1.2, '(\theta,\psi)','fontsize',14, 'FontWeight','normal','HorizontalAlignment','right','VerticalAlignment','top');

                app.Data.hs.quiverdazflat = quiver(0,0, 0*2 ,0*2,'color','k','linewidth',2);
                app.Data.hs.quiverdelflat = quiver(0,0, 0*2 ,0*2,'color','k','linewidth',2 );
                app.Data.hs.quivertvwflat = quiver(0,0, 0*2 ,0*2,'color',colors(3,:),'linewidth',2  ,'LineStyle',':');
                app.Data.hs.quivertvvflat = quiver(0,0, 0*2,0*2,'color',colors(4,:),'linewidth',2  ,'LineStyle',':');
                app.Data.hs.quivertvflat  = quiver(0,0, 0*2,0*2,'color','k','linewidth',2 ,'LineStyle','--');


                app.Data.hs.quivertJwflat = quiver(0,0, 0*2 ,0*2,'color',colors(3,:),'linewidth',1);
                app.Data.hs.quivertJvflat = quiver(0,0, 0*2 ,0*2,'color',colors(4,:),'linewidth',1);
                app.Data.hs.quivertAllvflat  = quiver(0,0, 0*2 ,0*2,'color','k','linewidth',1);

                legend([app.Data.hs.quivertJvflat app.Data.hs.quivertJwflat app.Data.hs.quivertAllvflat],{'Linear motion' 'Rotational motion' 'Total motion'})
            end

            % get data and update

            R = 1; % radius of the eye
            step = 10;
            [az, el] = meshgrid(deg2rad(-80:step:80),deg2rad(-80:step:80)); % azimuths and elevations to include

            paz = deg2rad(app.Values.Azimuth);
            pel = deg2rad(app.Values.Elevation);

            sys = app.Values.CoordinateSystem;
            set(app.Data.hs.TextSys, 'String',sys)

            % calculate spherical coordinates depending on the coordinate system
            switch(sys)
                case 'Fick'
                    [x,y,z] = Geometry3D.FickToSphere(az,el);
                    [px,py,pz] = Geometry3D.FickToSphere(paz,pel);

                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.FickLinearJacobian(az, el);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.FickRotationalJacobian(az, el);
                    Jv = cat(1, reshape([dazdx(:)'; dazdy(:)'; dazdz(:)'],1,3, length(az(:))),  reshape([deldx(:)'; deldy(:)' ;deldz(:)'],1,3, length(az(:))));
                    Jw = cat(1, reshape([rdazwx(:)'; rdazwy(:)'; rdazwz(:)'],1,3, length(az(:))),  reshape([rdelwx(:)' ;rdelwy(:)' ;rdelwz(:)'],1,3, length(az(:))));

                    [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.FickLinearInverseJacobian(paz, pel);
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.FickLinearJacobian(paz, pel);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.FickRotationalJacobian(paz, pel);

                case 'Helmholtz'
                    [x,y,z] = Geometry3D.HelmholtzToSphere(az,el);
                    [px,py,pz] = Geometry3D.HelmholtzToSphere(paz,pel);

                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HelmholtzLinearJacobian(az, el);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.HelmholtzRotationalJacobian(az, el);
                    Jv = cat(1, reshape([dazdx(:)'; dazdy(:)'; dazdz(:)'],1,3, length(az(:))),  reshape([deldx(:)'; deldy(:)' ;deldz(:)'],1,3, length(az(:))));
                    Jw = cat(1, reshape([rdazwx(:)'; rdazwy(:)'; rdazwz(:)'],1,3, length(az(:))),  reshape([rdelwx(:)' ;rdelwy(:)' ;rdelwz(:)'],1,3, length(az(:))));

                    
                    [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.HelmholtzLinearInverseJacobian(paz, pel);
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HelmholtzLinearJacobian(paz, pel);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.HelmholtzRotationalJacobian(paz, pel);
                case 'Harms'
                    [x,y,z] = Geometry3D.HarmsToSphere(az,el);
                    [px,py,pz] = Geometry3D.HarmsToSphere(paz,pel);

                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HarmsLinearJacobian(az, el);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.HarmsRotationalJacobian(az, el);
                    Jv = cat(1, reshape([dazdx(:)'; dazdy(:)'; dazdz(:)'],1,3, length(az(:))),  reshape([deldx(:)'; deldy(:)' ;deldz(:)'],1,3, length(az(:))));
                    Jw = cat(1, reshape([rdazwx(:)'; rdazwy(:)'; rdazwz(:)'],1,3, length(az(:))),  reshape([rdelwx(:)' ;rdelwy(:)' ;rdelwz(:)'],1,3, length(az(:))));

                    [dxdaz, dydaz, dzdaz, dxdel, dydel, dzdel] = Geometry3D.HarmsLinearInverseJacobian(paz, pel);
                    [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HarmsLinearJacobian(paz, pel);
                    [rdazwx, rdazwy, rdazwz, rdelwx, rdelwy, rdelwz] = Geometry3D.HarmsRotationalJacobian(paz, pel);
                case 'Hess'
                    [x,y,z] = Geometry3D.HessToSphere(az,el);
                    [px,py,pz] = Geometry3D.HessToSphere(paz,pel);
            end

            % scale by radius so the translations are in the right units (m)
            z = z*R;
            x = x*R;
            y = y*R;

            % udpate sphere
            set(app.Data.hs.meshSphere, 'xdata',x,'ydata',y,'zdata',z)


            w = deg2rad( [app.Values.AngularVelocityX_deg_s_ , app.Values.AngularVelocityY_deg_s_ , app.Values.AngularVelocityZ_deg_s_] );
            v = [app.Values.LinearVelocityX_m_s_, app.Values.LinearVelocityY_m_s_, app.Values.LinearVelocityZ_m_s_];
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
            tvw = Ji*[rdazwx, rdazwy, rdazwz; rdelwx, rdelwy, rdelwz]*w';
            tvv = Ji*[dazdx, dazdy, dazdz; deldx, deldy, deldz]*v';
            tv = tvw + tvv;

            set(app.Data.hs.quivertvw, 'xdata',px,'ydata',py,'zdata',pz);
            set(app.Data.hs.quivertvv, 'xdata',px,'ydata',py,'zdata',pz);
            set(app.Data.hs.quivertv, 'xdata',px,'ydata',py,'zdata',pz);
            set(app.Data.hs.quivertvw, 'UData',tvw(1),'VData',tvw(2),'WData',tvw(3));
            set(app.Data.hs.quivertvv, 'UData',tvv(1),'VData',tvv(2),'WData',tvv(3));
            set(app.Data.hs.quivertv, 'UData',tv(1),'VData',tv(2),'WData',tv(3));



            % update flat
%             set(app.Data.hs.meshflat, 'xdata',x*0,'ydata',y./x,'zdata',z./x)
%             set(app.Data.hs.meshpointflat, 'xdata',psx*0.1+0,'ydata',psy*0.05+py./px,'zdata',psz*0.05+pz./px)
%             set(app.Data.hs.textpointflat, 'Position', [0,py./px,pz./px])
%             set(app.Data.hs.meshflat, 'xdata',x*0,'ydata',az,'zdata',el)
%             set(app.Data.hs.meshpointflat, 'xdata',psx*0.1+0,'ydata',psy*0.05+paz,'zdata',psz*0.05+pel)
            set(app.Data.hs.textpointflat, 'Position', [paz,pel])


            set(app.Data.hs.quiverdazflat, 'xdata',paz,'ydata',pel);
            set(app.Data.hs.quiverdelflat, 'xdata',paz,'ydata',pel);
            set(app.Data.hs.quiverdazflat, 'UData',dydaz,'VData',dzdaz);
            set(app.Data.hs.quiverdelflat, 'UData',dydel,'VData',dzdel);

            tvw = [rdazwx, rdazwy, rdazwz; rdelwx, rdelwy, rdelwz]*w';
            tvv = [dazdx, dazdy, dazdz; deldx, deldy, deldz]*v';
            tv = tvw + tvv;


            ww = squeeze(pagemtimes(Jw,w'));
            vv = squeeze(pagemtimes(Jv,v'));
            allv = ww+vv;
            

            set(app.Data.hs.quivertJwflat, 'xdata', az(:),'ydata',el(:),     'UData',ww(1,:)' , 'VData',ww(2,:)');
            set(app.Data.hs.quivertJvflat, 'xdata', az(:), 'ydata',el(:),    'UData',vv(1,:)' , 'VData',vv(2,:)');
            set(app.Data.hs.quivertAllvflat, 'xdata', az(:),'ydata',el(:),   'UData',allv(1,:)' ,'VData',allv(2,:)');



            set(app.Data.hs.quivertvwflat, 'xdata',paz,'ydata',pel);
            set(app.Data.hs.quivertvvflat, 'xdata',paz,'ydata',pel);
            set(app.Data.hs.quivertvflat, 'xdata',paz,'ydata',pel);
            set(app.Data.hs.quivertvwflat, 'UData',tvw(1),'VData',tvw(2));
            set(app.Data.hs.quivertvvflat, 'UData',tvv(1),'VData',tvv(2));
            set(app.Data.hs.quivertvflat, 'UData',tv(1),'VData',tv(2));



            % draw flat sphere
            if ( 0)
                subplot(2,2,2,'nextplot','add');

                z = z*R;
                x = x*R;
                y = y*R;



                mesh(x*0,y,z,'FaceAlpha', 0.9,'facecolor',[1 1 1]);

                set(gca,'xtick',[],'ytick',[],'ztick',[])
                view(90,0)
                set(gca,'visible','off')


                line([R R ],[0 R*1.1 ],[0 0 ],'color',colors(1,:),'linewidth',2)
                line([R R ],[0 0 ],[0 R*1.1 ],'color',colors(5,:),'linewidth',2)
                hs = scatter3([R R ],[0 0 ],[0 0 ], 'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none' );
                text(1.1,0,0, 'x','FontWeight','normal','HorizontalAlignment','right','VerticalAlignment','top')
                text(0,1.4,0, 'y','FontWeight','normal','HorizontalAlignment','center')
                text(0,0,1.4, 'z','FontWeight','normal','HorizontalAlignment','center')

                [psx,psy,psz] = sphere(10);
                mesh(psx*0.05+px, psy*0.05+py,psz*0.05+pz, 'EdgeColor','none','FaceColor','k')
                text(px*1.2,py*1.2,pz*1.2, '(\theta,\psi)','fontsize',14, 'FontWeight','normal','HorizontalAlignment','left','VerticalAlignment','bottom')

            end
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
            el = asin( z / D );   % latitude
        end

        function [x, y, z] = HelmholtzToSphere(az,el)

            x = cos( az ) .* cos( el );
            y = sin( az );                 % latitude
            z = sin( el ) .* cos( az );   % longitude
        end

        function [az,el] = SphereToHelmholtz(x,y,z)
            D = sqrt( x.^2 + y.^2 + z.^2 );

            az = asin( y / D );  % latitude
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
            [dazdx, dazdy, dazdz, deldx, deldy, deldz] = Geometry3D.HelmholtzLinearJacobian(az, el);

            dxdaz = dazdx;
            dydaz = dazdy;
            dzdaz = dazdz;
            dxdel = deldx;
            dydel = deldy;
            dzdel = deldz;
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


        function [x, y, z] = HessToSphere(az, el)
            azdeg = rad2deg(az);

            z = sind(el);
            y = sind(az);
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

        function [az,el] = SphereToHess(~,y,z)
            D = sqrt( x.^2 + y.^2 + z.^2 );

            el = asin( z / D ); % latitude
            az = asin( y / D ); % latitude
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

    end
end

% TODO:
% https://work.thaslwanter.at/thLib/html/rotmat.html
