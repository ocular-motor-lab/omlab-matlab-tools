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

        function demoDisparity()
          
            app = InteractiveUI('Disparity Simulator',@(app) (Geometry3D.demoDisparityUpdate(app)), 0.2);
            app.AddDropDown('Stimulus',      1,  ["CROSS" "RANDOMPLANE" "GRIDPLANE" "HLINE" "VLINE"])
            app.AddSlider('IPD mm',           60, [10 100])
            app.AddSlider('Stimulus Scale',       1,  [0.1 10])
            app.AddSlider('Stimulus Distance',    40, [10 200])
            app.AddSlider('Fixation Distance',30, [10 200])
            app.AddSlider('Fixation X',        0, [-100 100])
            app.AddSlider('Fixation Y',        0, [-100 100])
            app.AddSlider('Torsion Version',  0,  [-20 20])
            app.AddSlider('Torsion Vergence', 0,  [-20 20])
            app.AddSlider('Plane slant',     0,  [-90 90])
            app.AddSlider('Plane Tilt',      0,  [0 90])
            app.AddDropDown('View3D',      1,  ["Oblique" "TOP" "SIDE"])


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


        function demoMotionFlow()
          
            app = InteractiveUI('Motion Flow Simulator',@(app) (Geometry3D.demoDisparityUpdate(app)), 0.2);
            app.AddDropDown('Stimulus',      1,  ["Ground plane" "Cloud"])
            app.AddSlider('Eye height',           100,[0    500])
            % app.AddSlider('Stimulus Scale',       1,  [0.1  10])
            app.AddSlider('Fixation Distance',    30, [10   200])
            app.AddSlider('Fixation X',           0,  [-100 100])
            app.AddSlider('Fixation Y',           0,  [-100 100])
            % app.AddSlider('Torsion Version',      0,  [-20  20])
            % app.AddSlider('Torsion Vergence',     0,  [-20  20])
            app.AddSlider('Ground plane slant',          0,  [-90  90])
            app.AddSlider('Ground plane tilt',           0,  [0    90])
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
            app = InteractiveUI('Quaternion demo',@(app) (Geometry3D.demoQuaternionUpdate(app)), 0.2);
            app.AddSlider('q0',1, [-1 1])
            app.AddSlider('q1',0, [-1 1])
            app.AddSlider('q2',0, [-1 1])
            app.AddSlider('q3',0, [-1 1])

            % app.AddDropDown('2D HVT Coordinate system',      1,  ["Helmholtz" "Fick" "Harns" "Hess"])
            app.AddDropDown('2D HVT Coordinate system',      1,  ["Helmholtz" "Fick"])
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

            % Position of eyes
            eyes.R.X = ipdCm/2;
            eyes.R.Y = 0;
            eyes.R.Z = 0;
            eyes.L.X = - ipdCm/2;%
            eyes.L.Y = 0;
            eyes.L.Z = 0;

            % Orientation of eyes in helmoltz coordinates, they make more
            % sense for 
            % right is positive, Up is positive
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


            t =table();

            % The eye coordinate system is right handed
            % with x pointing out
            % y pointing to left
            % z pointing up
            %
            % it is a bit messy to have a different coordinate system than 
            % for the 3D scene, but I can't think of the rotations otherway
            %

            % Get the new xs in cm for the two eyes (transform the dots from head reference
            % to eye reference [nothing about eye angle here!])
            t.RY = -(points.X - eyes.R.X);
            t.RZ = -points.Y;
            t.RX = points.Z;
            t.LY = -(points.X - eyes.L.X);
            t.LZ = -points.Y;
            t.LX = points.Z;

            % Rotate the points to be in eye reference frame
            t{:,{'RX' 'RY' 'RZ'  }} =  t{:,{'RX' 'RY' 'RZ'  }}*eyes.R.RM;
            t{:,{'LX' 'LY' 'LZ'  }} =  t{:,{'LX' 'LY' 'LZ'  }}*eyes.L.RM;

            % helmholtz coordinates make more sense for disparity
            % because the vertical rotation first does not change the
            % plane of regard
            t.RV = atan2d(t.RZ, t.RX);
            t.RH = atan2d(-t.RY, t.RX./abs(cosd(t.RV)));
            t.LV = atan2d(t.LZ, t.LX);
            t.LH = atan2d(-t.LY, t.LX./abs(cosd(t.LV)));
            badidx = abs(t.RV)>=90 | abs(t.RH)>=90 ;
            t.RV(badidx) = nan;
            t.RH(badidx) = nan;
            badidx = abs(t.LV)>=90 | abs(t.LH)>=90 ;
            t.LV(badidx) = nan;
            t.LH(badidx) = nan;

            % Get disparity!! % TODO this should be 3D disparity
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

        function [f, heyes, hfix, hscreen, hpoints, hspoints, hLRpoints, hdisparity, hLscreen] = demoDisparityInitPlots(points, eyePoints, screenPoints, eyes, leftEyeScreen, rightScreen)

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
                        [X,Y] = meshgrid([-(sizeStimCm/2):2:(sizeStimCm/2)],[-(sizeStimCm/2):2:(sizeStimCm/2)]);
                        worldPoints = table();
                        worldPoints.X = X(:);
                        worldPoints.Y = Y(:);
                        worldPoints.Z = zeros(size(worldPoints.X));

                    case 'HLINE'
                        % TODO: need to update the update funciton so the Z is not overwritten by a
                        % constant distance. Just shift things.
                        [X,Y] = meshgrid([-(sizeStimCm/2):0.5:(sizeStimCm/2)],[0]*ones(size([-(sizeStimCm/2):0.5:(sizeStimCm/2)])));
                        worldPoints = table();
                        worldPoints.X = X(:);
                        worldPoints.Y = Y(:);
                        worldPoints.Z = zeros(size(worldPoints.X));
                    case 'VLINE'
                        % TODO: need to update the update funciton so the Z is not overwritten by a
                        % constant distance. Just shift things.
                        [X,Y] = meshgrid([0]*ones(size([-(sizeStimCm/2):0.5:(sizeStimCm/2)])),[-(sizeStimCm/2):0.5:(sizeStimCm/2)]);
                        worldPoints = table();
                        worldPoints.X = X(:);
                        worldPoints.Y = Y(:);
                        worldPoints.Z = zeros(size(worldPoints.X));



                    case 'CROSS'
                        % TODO: need to update the update funciton so the Z is not overwritten by a
                        % constant distance. Just shift things.
                        [X,Y] = meshgrid([-(sizeStimCm/2):0.5:(sizeStimCm/2)],[0]*ones(size([-(sizeStimCm/2):0.5:(sizeStimCm/2)])));
                        worldPoints = table();
                        worldPoints.X = X(:);
                        worldPoints.Y = Y(:);
                        worldPoints.Z = zeros(size(worldPoints.X));
                        % TODO: need to update the update funciton so the Z is not overwritten by a
                        % constant distance. Just shift things.
                        [X,Y] = meshgrid([0]*ones(size([-(sizeStimCm/2):0.5:(sizeStimCm/2)])),[-(sizeStimCm/2):0.5:(sizeStimCm/2)]);
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
            h.sphere = surf(R*X,R*Y,R*Z,'FaceAlpha', alpha,'facecolor',[0.8 0.8 0.8],'edgecolor',0.2*[0.8 0.8 0.8])
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

            [xcirclesX, xcirclesY] = meshgrid(steps,points);
            [ycirclesX, ycirclesY] = meshgrid(points,steps);

            
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



            switch(app.Values.x2DHVTCoordinateSystem)
                case 'Fick'
                    z = sind(xcirclesX);
                    y = sind(xcirclesY).*cosd(xcirclesX);
                    x = cosd(xcirclesX).*cosd(xcirclesY);
                case 'Helmholtz'
                    y = sind(xcirclesX);
                    x = sind(xcirclesY).*cosd(xcirclesX);
                    z = cosd(xcirclesX).*cosd(xcirclesY);
                case 'Listings(AET)'
                    x = sind(xcirclesX);
                    z = sind(xcirclesY).*cosd(xcirclesX);
                    y = cosd(xcirclesX).*cosd(xcirclesY);
                    
                % case 'Harns'
                %     y = sind(xcirclesX).*cosd(xcirclesY);
                %     x = sind(xcirclesY);
                %     z = cosd(xcirclesX).*cosd(xcirclesY);
                % case 'Hess'
                %     y = sind(xcirclesX);
                %     x = sind(xcirclesY);
                %     z = cosd(xcirclesX).*cosd(xcirclesY);
            end

            set(app.Data.h.sphere, 'xdata',x, 'ydata',y, 'zdata',z)

            
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
end

% TODO:
% https://work.thaslwanter.at/thLib/html/rotmat.html