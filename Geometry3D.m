classdef Geometry3D
    %Geometry3D Summary of this class goes here
    %   Detailed explanation goes here
    %
    %

    properties
    end

    methods(Static)

        function demoDisparity(stimType)
            if (~exist("stimType", "var"))
                stimType = 'GRID';
            end
            %%
            sizeStimCm = 20;


            % Make some dots in world coordinates (XYZ)
            % 
            % Make a table to start putting everything together
            switch (stimType)
                case 'PLANE'
                    numDots = 100;
                    worldPoints = table();
                    worldPoints.X = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
                    worldPoints.Y = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
                    worldPoints.Z = zeros(numDots, 1);
                case 'GRID'
                    % TODO: need to update the update funciton so the Z is not overwritten by a
                    % constant distance. Just shift things.
                    [X,Y] = meshgrid([-10:2:10],[-10:2:10]);
                     worldPoints = table();
                     worldPoints.X = X(:);
                     worldPoints.Y = Y(:);
                     worldPoints.Z = zeros(size(worldPoints.X));

                case 'HLINE'
                    % TODO: need to update the update funciton so the Z is not overwritten by a
                    % constant distance. Just shift things.
                    [X,Y] = meshgrid([-10:0.5:10],[0]*ones(size([-10:0.5:10])));
                     worldPoints = table();
                     worldPoints.X = X(:);
                     worldPoints.Y = Y(:);
                     worldPoints.Z = zeros(size(worldPoints.X));
                case 'VLINE'
                    % TODO: need to update the update funciton so the Z is not overwritten by a
                    % constant distance. Just shift things.
                    [X,Y] = meshgrid([0]*ones(size([-10:0.5:10])),[-10:0.5:10]);
                     worldPoints = table();
                     worldPoints.X = X(:);
                     worldPoints.Y = Y(:);
                     worldPoints.Z = zeros(size(worldPoints.X));



                case 'CROSS'
                    % TODO: need to update the update funciton so the Z is not overwritten by a
                    % constant distance. Just shift things.
                    [X,Y] = meshgrid([-10:0.5:10],[0]*ones(size([-10:0.5:10])));
                     worldPoints = table();
                     worldPoints.X = X(:);
                     worldPoints.Y = Y(:);
                     worldPoints.Z = zeros(size(worldPoints.X));
                    % TODO: need to update the update funciton so the Z is not overwritten by a
                    % constant distance. Just shift things.
                    [X,Y] = meshgrid([0]*ones(size([-10:0.5:10])),[-10:0.5:10]);
                     worldPoints2 = table();
                     worldPoints2.X = X(:);
                     worldPoints2.Y = Y(:);
                     worldPoints2.Z = zeros(size(worldPoints2.X));

                     worldPoints = cat(1,worldPoints, worldPoints2);
            end

            app = InteractiveUI('DisparityCalculator',@(app) (Geometry3D.demoDisparityUpdate(app, worldPoints)), 0.2);
            app.AddControl('fixationDistance',30, [10 200])
            app.AddControl('stimDistance',    40, [10 200])
            app.AddControl('stimScale',       1,  [0.1 10])
            app.AddControl('ipdmm',           60, [10 100])
            app.AddControl('torsionVersion',  0,  [-20 20])
            app.AddControl('torsionVergence', 0,  [-20 20])
            app.AddControl('Plane slant',     0,  [-90 90])
            app.AddControl('Plane Tilt',      0,  [0 90])


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

        function demoQuaternion()
            app = InteractiveUI('Quaternion demo',@(app) (Geometry3D.demoQuaternionUpdate(app)), 0.2);
            app.AddControl('q0',1, [-1 1])
            app.AddControl('q1',0, [-1 1])
            app.AddControl('q2',0, [-1 1])
            app.AddControl('q3',0, [-1 1])

            app.AddControl('H',0, [-90 90])
            app.AddControl('V',0, [-90 90])
            app.AddControl('T',0, [-90 90])

            app.AddControl('x',0, [-1 1])
            app.AddControl('y',0, [-1 1])
            app.AddControl('z',0, [-1 1])
            app.AddControl('a',0, [-180 180])
            

            app.Open();
        end

        function eyes = MakeEyes(ipdCm, fixationSpot, torsionVersion, torsionVergence)

            % Position of eyes
            eyes.R.X =  ipdCm/2;
            eyes.R.Y = 0;
            eyes.R.Z = 0;
            eyes.L.X = -ipdCm/2;%
            eyes.L.Y = 0;
            eyes.L.Y = 0;

            % Orientation of eyes % TODO: DEFINE COORDINATE SYSTEM right
            % now , is it fick? 
            eyes.R.H = -atan2d(eyes.R.X, fixationSpot.Z);
            eyes.R.V = 0;
            eyes.R.T = torsionVersion - torsionVergence;
            eyes.L.H = -atan2d(eyes.L.X, fixationSpot.Z);
            eyes.L.V = 0;
            eyes.L.T = torsionVersion + torsionVergence;


            % Get the rotation matrices describing the eye angles
            REyeRMx = Geometry3D.RotZ(deg2rad(eyes.R.H)); % all rows are the same so we can just use the first one
            REyeRMy = Geometry3D.RotY(deg2rad(eyes.R.V));
            REyeRMz = Geometry3D.RotX(deg2rad(eyes.R.T));
            eyes.R.RM = REyeRMx*REyeRMy*REyeRMz; % Rot matrix defining the orientation of the eye relative to straight ahead direction.

            LEyeRMx = Geometry3D.RotZ(deg2rad(eyes.L.H)); % all rows are the same so we can just use the first one
            LEyeRMy = Geometry3D.RotY(deg2rad(eyes.L.V));
            LEyeRMz = Geometry3D.RotX(deg2rad(eyes.L.T));
            eyes.L.RM = LEyeRMx*LEyeRMy*LEyeRMz; % Rot matrix defining the orientation of the eye relative to straight ahead direction.

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


            t =table();


            % Get the new xs in cm for the two eyes (transform the dots from head reference
            % to eye reference [nothing about eye angle here!])
            t.RX = points.X - eyes.R.X; % both in cm
            t.RY = points.Y;
            t.RZ = points.Z;
            t.LX = points.X - eyes.L.X; % no need to do this for y bc the two eyes are only displaced horizontally on the head
            t.LY = points.Y;
            t.LZ = points.Z;

            % Get the ANGLES of the dots relative to the (straight ahead) eye position in cm
            % TODO: do proper quaternions/rot matrices here maybe
            t.DotsAngleRx = atan2d(t.RX, points.Z);
            t.DotsAngleLx = atan2d(t.LX, points.Z);
            t.DotsAngleRy = atan2d(points.Y, (points.Z)./abs(cosd(t.DotsAngleRx)));
            t.DotsAngleLy = atan2d(points.Y, (points.Z)./abs(cosd(t.DotsAngleLx)));

            % % % % % % Calculate the angles of the dots in the eye reference frame accounting
            % % % % % % for the angles (vergence) of the eyes
            % % % % % % TODO: do proper quaternions/rot matrices here maybe
            % % % % % t.RAnglex = t.DotsAngleRx - eyes.R.H; % both in deg
            % % % % % t.LAnglex = t.DotsAngleLx - eyes.L.H;
            % % % % % t.RAngley = t.DotsAngleRy;
            % % % % % t.LAngley = t.DotsAngleLy;


            % Get the rotation matrices describing the eye angles (X coming
            % out of the eye)
            REyeRM = eyes.R.RM;
            LEyeRM = eyes.L.RM;

            % Calculate the angles of the dots in the eye reference frame accounting
            % for the angles (vergence) of the eyes
            % for i = 1:height(t)
            %     RDotRMx = Geometry3D.RotZ(deg2rad(t.DotsAngleRx(i)));
            %     RDotRMy = Geometry3D.RotY(deg2rad(t.DotsAngleRy(i)));
            %     t.RDotRM{i} = RDotRMx*RDotRMy; % Rot matrix defining the orientation (angles) of a dot relative to straight ahead eye direction.
            %     t.RDotinEyeRM{i} = (t.RDotRM{i}'*REyeRM)'; % rot matrix defining the orientation of the dot relative to the eye actual direction.
            % 
            %     LDotRMx = Geometry3D.RotZ(deg2rad(t.DotsAngleLx(i)));
            %     LDotRMy = Geometry3D.RotY(deg2rad(t.DotsAngleLy(i)));
            %     t.LDotRM{i} = LDotRMx*LDotRMy; % Rot matrix defining the orientation (angles) of a dot relative to straight ahead eye direction.
            %     t.LDotinEyeRM{i} = (t.LDotRM{i}'*LEyeRM)'; % rot matrix defining the orientation of the dot relative to the eye actual direction.
            % end

            % TODO: Document!!! and make sure X and Y are correctly
            % flipped. I Don't quite understand why. 
            % Rotate the points to be in eye reference frame
            Rzyx =  t{:,{'RZ' 'RX' 'RY'  }}*REyeRM;
            Lzyx =  t{:,{'LZ' 'LX' 'LY'  }}*LEyeRM;

            % % Get the hor vert tor rotations in fick
            % for i = 1:height(t)
            %     RHVT = Geometry3D.RotMat2Fick(t.RDotinEyeRM{i});
            %     LHVT = Geometry3D.RotMat2Fick(t.LDotinEyeRM{i});
            % 
            %     t.RH(i) = rad2deg(RHVT(1));
            %     t.RV(i) = rad2deg(RHVT(2));
            %     t.LH(i) = rad2deg(LHVT(1));
            %     t.LV(i) = rad2deg(LHVT(2));
            % end

            t.RH = atan2d(Rzyx(:,2), Rzyx(:,1));
            t.RV = atan2d(Rzyx(:,3), (Rzyx(:,1))./abs(cosd(Rzyx(:,2))));
            t.LH = atan2d(Lzyx(:,2), Lzyx(:,1));
            t.LV = atan2d(Lzyx(:,3), (Lzyx(:,1))./abs(cosd(Lzyx(:,2))));

            % Get disparity!!
            t.HDisparity = t.LH - t.RH;
            t.VDisparity = t.LV - t.RV;

            eyepoints = t;
        end


        function screenPoints  = PointsEyesToScreen(eyes, eyePoints, LeyeScreen, ReyeScreen)

            % TODO: correct for fick/helm or whatever coordinate systems we
            % have in the eye posints
            
            % TODO: Review how eye position is added (probably also needs
            % to be a rotation

            % This needs to be different if we are simulating torsion or
            % not

            screenPoints.LX = LeyeScreen.middleX + LeyeScreen.ScreenDistance*tand(eyes.L.H + eyePoints.LH)*LeyeScreen.pixPerCmWidth;
            screenPoints.LY = LeyeScreen.middleY + LeyeScreen.ScreenDistance*tand(eyes.L.V + eyePoints.LV)*LeyeScreen.pixPerCmHeight;
            screenPoints.RX = ReyeScreen.middleX + ReyeScreen.ScreenDistance*tand(eyes.R.H + eyePoints.RH)*ReyeScreen.pixPerCmWidth;
            screenPoints.RY = ReyeScreen.middleY + ReyeScreen.ScreenDistance*tand(eyes.R.V + eyePoints.RV)*ReyeScreen.pixPerCmHeight;
                
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

        function [f, heyes, hfix, hpoints, hspoints, hLRpoints, hdisparity, hLscreen] = demoDisparityInit(points, eyePoints, screenPoints, eyes)

            t = points(1:end-1,:);
            fp = points(end,:);

            % View the dots, eyes, and fixation dot
            scr_siz = get(0,'ScreenSize');
            margin = floor(0.1*(scr_siz(4)));
            f = figure('position',floor([...
                margin...
                margin...
                scr_siz(3)*2.8/4 ...
                scr_siz(4)*2/4 ...
                ]));
            subplot(2,4,[1 2 5 6], 'nextplot','add')

            hpoints=plot3(t.X,t.Z, t.Y,'o','Color','k');
            hspoints=plot3(t.X,t.Z, zeros(size(t.Y)),'o','Color',[0.8 0.8 0.8]);
            hfix=plot3(fp.X,fp.Z,fp.Y,'o','Color','r','LineWidth',2, 'markersize',20); % fixation spot
            heyes.c = plot3([eyes.L.X eyes.R.X], [eyes.L.Y eyes.R.Y], [0 0],'o','Color','b', 'markersize',50); % left eye fixation spot and right eye
            heyes.l = plot3([2*[-1 0 1 0  0 0 0] eyes.L.X 0 eyes.R.Y], [2*[ 0 0 0 0 -1 0 1] eyes.L.Y 0 eyes.R.Y],zeros(1,10),'-','Color','b'); % left eye fixation spot and right eye

            xlim(1.2*[min([t.X;t.Y]) max([t.X;t.Y])])
            zlim(1.2*[min([t.X;t.Y]) max([t.X;t.Y])])
            % ylim([-5 min(100,max(fp.Z, max(t.Z)))*2])
            ylim([-5 200])
            xlabel('X (cm)')
            zlabel('Y (cm)')
            ylabel('Z (cm)')
            view(-75,10)
            title('3D world');

            hLRpoints = [];
            t = eyePoints;
            subplot(2,4,3)
            % hLRpoints.R = plot(t.RH, t.RV,'o');
            ra = atan2(t.RV, t.RH );
            re = sqrt(t.RH.^2 + t.RV.^2);
            hLRpoints.R = polarplot(ra, re,'o');
            set(gca, 'nextplot','add')

            ra = atan2(t.LV, t.LH );
            re = sqrt(t.LH.^2 + t.LV.^2);
            hLRpoints.L = polarplot(ra, re,'o');
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
            title('Retinal position');

            subplot(2,4,4)
            hdisparity = quiver(t.RH, t.RV, t.LH-t.RH, t.LV-t.RV, 'AutoScale', "off");
            % for i = 1:height(t)
            %     line([t.RAnglex(i,1) t.LAnglex(i,1)],[t.RAngley(i,1) t.LAngley(i,1)])
            % end
            grid
            title('Disparity')
            set(gca,'xlim',[-30 30],'ylim',[-30 30])

            % subplot(2,4,7)
            % for i = 1:height(t)
            %     line([t.RH(i,1) t.LH(i,1)],[t.RV(i,1) t.LV(i,1)])
            % end
            % title('Disparity With Torsion')

            subplot(2,4,7, 'nextplot','add')
            hLscreen.LPoints = plot(screenPoints.LX(1:end-1), screenPoints.LY(1:end-1), '+');
            
            hLscreen.LFP = plot(screenPoints.LX(end), screenPoints.LY(end), 'ro');
            grid
            set(gca,'xlim',[0 1920],'ylim',[0 1080])
            set(gca,'PlotBoxAspectRatio',[16 9 1])
            title('Left eye screen (sim torsion)');


            subplot(2,4,8, 'nextplot','add')
            hLscreen.RPoints = plot(screenPoints.RX(1:end-1), screenPoints.RY(1:end-1), '+');
            
            hLscreen.RFP = plot(screenPoints.RX(end), screenPoints.RY(end), 'ro');
            grid
            set(gca,'xlim',[0 1920],'ylim',[0 1080])
            set(gca,'PlotBoxAspectRatio',[16 9 1])
            title('Right eye screen (sim torsion)');
            
        end

        function demoDisparityUpdate(app, worldPoints)

            sliderValues = app.SliderValues;

            %% Update the points, eyes and screen acordint to the sliders
            app.Data.FixationSpot.Z = sliderValues.fixationDistance;

            leftEyeScreen = Geometry3D.MakeScreen(app.Data.Screen.SizeCm, app.Data.Screen.ResPix, app.Data.Screen.Distance, app.Data.Screen.Slant);
            rightEyeScreen = Geometry3D.MakeScreen(app.Data.Screen.SizeCm, app.Data.Screen.ResPix, app.Data.Screen.Distance, app.Data.Screen.Slant);
            
            eyes = Geometry3D.MakeEyes(sliderValues.ipdmm/10, app.Data.FixationSpot, sliderValues.torsionVersion, sliderValues.torsionVergence);

            % Rotate the world points to apply the slant (rotates around
            % 0,0,0)
            worldPoints.X = worldPoints.X*sliderValues.stimScale;
            worldPoints.Y = worldPoints.Y*sliderValues.stimScale;
            worldPoints.Z = zeros(size(worldPoints.Z));
            
            % Rotate by slant and tilt
            R = Geometry3D.Quat2RotMat(Geometry3D.AxisAngle2Quat([cosd(sliderValues.PlaneTilt) sind(sliderValues.PlaneTilt) 0],deg2rad(sliderValues.PlaneSlant)));
            worldPoints{:,:} = (R*worldPoints{:,:}')';

            % Displace by distance
            worldPoints.Z = worldPoints.Z + sliderValues.stimDistance;

            % Add the fixation spot to the world points to conver to eye
            % and screen points.
            worldPoints{end+1,:} = [app.Data.FixationSpot.X app.Data.FixationSpot.Y app.Data.FixationSpot.Z];

            %% Get the points in eye and screen coordinates
            eyePoints = Geometry3D.Points3DToEyes(worldPoints, eyes);
            screenPoints = Geometry3D.PointsEyesToScreen(eyes, eyePoints, leftEyeScreen, rightEyeScreen);


            if ( ~isfield(app.Data, "f") || ~isvalid(app.Data.f))
                % If figure does not exist create it with all the plots and
                % handles to them
                [f, heyes, hfix, hpoints, hspoints, hLRpoints, hdisparity, hscreen] = Geometry3D.demoDisparityInit(worldPoints, eyePoints, screenPoints, eyes);
                app.Data.f = f;
                app.Data.heyes = heyes;
                app.Data.hfix = hfix;
                app.Data.hpoints = hpoints;
                app.Data.hspoints = hspoints;
                app.Data.hLRpoints = hLRpoints;
                app.Data.hdisparity = hdisparity;
                app.Data.hscreen = hscreen;
            end
        
            % update 3D plot
            set(app.Data.hfix, 'ydata', sliderValues.fixationDistance);
            set(app.Data.heyes.l, 'xdata', [eyes.L.X 0 eyes.R.X], 'ydata', [0 sliderValues.fixationDistance 0], 'zdata', [0 0 0]);
            set(app.Data.heyes.c, 'xdata', [eyes.L.X eyes.R.X], 'ydata', [0 0], 'zdata', [0 0]);
            set(app.Data.hpoints, 'xdata', worldPoints.X, 'ydata',worldPoints.Z, 'zdata', worldPoints.Y);
            set(app.Data.hspoints, 'xdata', worldPoints.X, 'ydata',worldPoints.Z, 'zdata', app.Data.hspoints.Parent.XLim(1)*ones(size(worldPoints.Y)));


            % update retina plot
            ra = atan2(eyePoints.RV, eyePoints.RH );
            re = sqrt(eyePoints.RH.^2 + eyePoints.RV.^2);
            set(app.Data.hLRpoints.L, 'ThetaData', ra, 'RData', re);
            ra = atan2(eyePoints.LV, eyePoints.LH );
            re = sqrt(eyePoints.LH.^2 + eyePoints.LV.^2);
            set(app.Data.hLRpoints.R, 'ThetaData', ra, 'RData', re);
            
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
        function [f, h] = demoQuaternionInit(app)

            f = figure;

            h = axis;

        end

        function demoQuaternionUpdate(app)
            
            sliderValues = app.SliderValues;
            q = [sliderValues.q0 sliderValues.q1 sliderValues.q2 sliderValues.q3];

            if ( ~isfield(app.Data, "f") || ~isvalid(app.Data.f))
                % If figure does not exist create it with all the plots and
                % handles to them
                [f, h] = Geometry3D.demoQuaternionInit();
                app.Data.f = f;
                app.Data.h = h;

                app.Data.LastQ = q;
            end

            
            q1 = app.Data.LastQ;

            % Update only the 3 values not changed by the user
            % so as to keep a valid quaternion. If the other three values
            % are zero transform to aligned with 1 0 0 because otherwise
            % there is a singularity
            c = (q-q1) == 0;
            if ( sum(c)==3)
                nc = q(c);
                if ( sum(nc) == 0)
                    nc = [1 0 0];
                end

                q(c) = nc*sqrt(1-q(~c)^2)/norm(nc);

                app.Data.LastQ = q;
                app.SliderValues.q0  = q(1);
                app.SliderValues.q1  = q(2);
                app.SliderValues.q2  = q(3);
                app.SliderValues.q3  = q(4);
            end

            [axis,angle] = Geometry3D.Quat2AxisAngle(q);
            HVT = Geometry3D.RotMat2Fick( Geometry3D.Quat2RotMat(q));

            app.SliderValues.H  = rad2deg(HVT(1));
            app.SliderValues.V  = rad2deg(HVT(2));
            app.SliderValues.T  = rad2deg(HVT(3));
            app.SliderValues.x  = axis(1);
            app.SliderValues.y  = axis(2);
            app.SliderValues.z  = axis(3);
            app.SliderValues.a  = angle;

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

        % HORIZONTAL ROTATION
        function M = RotZ(theta)
            M = [   cos(theta)  -sin(theta)     0;
                sin(theta)  cos(theta)      0;
                0           0               1];
        end

        % VERTICAL ROTATION
        function M = RotY(phi)
            M = [   cos(phi)  0               sin(phi);
                0           1               0;
                -sin(phi) 0               cos(phi)];
        end

        % TORSIONAL ROTATION
        function M = RotX(psi)
            M = [   1           0               0;
                0           cos(psi)      -sin(psi);
                0           sin(psi)      cos(psi)];
        end
    end

    methods(Static)

        % Fick to rotation matrix
        function M = Fick2RotMat(HVT)

            theta = HVT(1);
            phi = HVT(2);
            psi = HVT(3);

            M = Geometry3D.RotZ(theta)*Geometry3D.RotY(phi)*Geometry3D.RotX(psi);
        end

        % Rotation matrix to Fick
        function HVT = RotMat2Fick(M)
            % TODO: Double check this.
            Rzx = M(3,1);
            Ryx = M(2,1);
            Rzy = M(3,2);

            HVT(2) = -asin(Rzx);
            HVT(1) = asin(Ryx/cos(HVT(2)));
            HVT(3) = asin(Rzy/cos(HVT(2)));
        end

        function M = Helm2RotMat(HVT)
            theta = HVT(1);
            phi = HVT(2);
            psi = HVT(3);

            M = Geometry3D.RotY(phi)*Geometry3D.RotZ(theta)*Geometry3D.RotX(psi);
        end

        function HVT = RotMat2Helm(HVT)
            Rzx = M(3,1);
            Ryx = M(2,1);
            Rzy = M(3,2);

            HVT(1) = -asin(Rzx);
            HVT(2) = asin(Ryx/cos(HVT(1)));
            HVT(3) = asin(Rzy/cos(HVT(1)));
        end

        function M = List2Mat(ADT)
            A = ADT(1);
            D = ADT(2);
            T = ADT(3);
            M = Geometry3D.RotX(A)*Geometry3D.RotZ(D)*Geometry3D.RotX(T-A);
        end

        function q = RotMat2Quat(M)
            % From quaternion navy book
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
end

% TODO:
% https://work.thaslwanter.at/thLib/html/rotmat.html