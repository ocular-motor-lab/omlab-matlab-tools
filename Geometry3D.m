classdef Geometry3D
    %Geometry3D Summary of this class goes here
    %   Detailed explanation goes here

    properties
    end

    methods(Static)

        function demoDisparity(stimType)
            if (~exist("stimType"))
                stimType = 'FLAT';
            end
            %%
            numDots = 90;
            sizeStimCm = 20;
            stimDistance = 10;


            % Make some dots in world coordinates (XYZ)
            % 
            % Make a table to start putting everything together
            switch (stimType)
                case 'FLAT'
                    worldPoints = table();
                    worldPoints.X = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
                    worldPoints.Y = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
                    worldPoints.Z = ones(numDots, 1)*stimDistance;
                case 'SLANT'
                    % TODO: need to update the update funciton so the Z is not overwritten by a
                    % constant distance. Just shift things.
                    worldPoints = table();
                    worldPoints.X = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
                    worldPoints.Y = rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
                    worldPoints.Z = ones(numDots, 1)*stimDistance + worldPoints.X;
            end


            idpCm = 6;
            fixationDistance = 400;
            torsionVersion = 30;
            torsionVergence = 0;


            fixationSpot = [];
            fixationSpot.X = 0;
            fixationSpot.Y = 0;  
            fixationSpot.Z = fixationDistance;

            worldPoints{end+1,:} = [fixationSpot.X fixationSpot.Y fixationSpot.Z];

            eyes = Geometry3D.MakeEyes(idpCm, fixationSpot, torsionVersion, torsionVergence);

            sizeCm = [30*16/9 30];
            resPix = [1920 1080];
            distance = 57;
            slant = 0;
            LeyeScreen = Geometry3D.MakeScreen(sizeCm, resPix, distance, slant);
            ReyeScreen = Geometry3D.MakeScreen(sizeCm, resPix, distance, slant);

            eyePoints = Geometry3D.Points3DToEyes(worldPoints, eyes);

            screenPoints = Geometry3D.PointsEyesToScreen(eyes, eyePoints, LeyeScreen, ReyeScreen);

            [heyes, hfix, hpoints, hLRpoints, hdisparity, hscreen] = Geometry3D.PlotPoints(worldPoints, eyePoints, screenPoints, eyes, fixationDistance);

            s = InteractiveUI('DisparityCalculator',@(sliderValues) (Geometry3D.demoDisparityUpdate(sliderValues, heyes, hfix, hpoints, hLRpoints, hdisparity, hscreen, worldPoints)), 0.2);
            s.AddControl('fixationDistance',30, [10 200])
            s.AddControl('stimDistance',    40, [10 200])
            s.AddControl('ipdmm',           60, [10 100])
            s.AddControl('torsionVersion',  0,  [0 30])
            s.AddControl('torsionVergence', 0,  [0 10])
            s.AddControl('nDots',           90, [0 100])
            s.AddControl('sizeStimCm',      10, [0 100])
            s.Open();

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

            screen.distance = distance;
            screen.distance = distance;
            
            screen.slant = slant;
        end


        function [eyepoints eyes] = Points3DToEyes(points, eyes )


            t =table();


            % Get the new xs in cm for the two eyes (transform the dots from head reference
            % to eye reference [nothing about eye angle here!])
            t.newxR = points.X - eyes.R.X; % both in cm
            t.newxL = points.X - eyes.L.X; % no need to do this for y bc the two eyes are only displaced horizontally on the head

            % Get the ANGLES of the dots relative to the (straight ahead) eye position in cm
            % TODO: do proper quaternions/rot matrices here maybe
            t.DotsAngleRx = atan2d(t.newxR, points.Z);
            t.DotsAngleLx = atan2d(t.newxL, points.Z);
            t.DotsAngleRy = atan2d(points.Y, (points.Z)./abs(cosd(t.DotsAngleRx)));
            t.DotsAngleLy = atan2d(points.Y, (points.Z)./abs(cosd(t.DotsAngleLx)));

            % % % % % % Calculate the angles of the dots in the eye reference frame accounting
            % % % % % % for the angles (vergence) of the eyes
            % % % % % % TODO: do proper quaternions/rot matrices here maybe
            % % % % % t.RAnglex = t.DotsAngleRx - eyes.R.H; % both in deg
            % % % % % t.LAnglex = t.DotsAngleLx - eyes.L.H;
            % % % % % t.RAngley = t.DotsAngleRy;
            % % % % % t.LAngley = t.DotsAngleLy;


            % Get the rotation matrices describing the eye angles
            REyeRM = eyes.R.RM;
            LEyeRM = eyes.L.RM;

            % Calculate the angles of the dots in the eye reference frame accounting
            % for the angles (vergence) of the eyes
            for i = 1:height(t)
                RDotRMx = Geometry3D.RotZ(deg2rad(t.DotsAngleRx(i)));
                RDotRMy = Geometry3D.RotY(deg2rad(t.DotsAngleRy(i)));
                t.RDotRM{i} = RDotRMx*RDotRMy; % Rot matrix defining the orientation (angles) of a dot relative to straight ahead eye direction.
                t.RDotinEyeRM{i} = (t.RDotRM{i}'*REyeRM)'; % rot matrix defining the orientation of the dot relative to the eye actual direction.

                LDotRMx = Geometry3D.RotZ(deg2rad(t.DotsAngleLx(i)));
                LDotRMy = Geometry3D.RotY(deg2rad(t.DotsAngleLy(i)));
                t.LDotRM{i} = LDotRMx*LDotRMy; % Rot matrix defining the orientation (angles) of a dot relative to straight ahead eye direction.
                t.LDotinEyeRM{i} = (t.LDotRM{i}'*LEyeRM)'; % rot matrix defining the orientation of the dot relative to the eye actual direction.
            end

            % test = pagetranspose(pagemtimes(pagetranspose(cat(3,t.RDotRM{:})),REyeRM));
            %
            % test = pagemtimes(t.RDotRM,REyeRM)

            % Get the hor vert tor rotations in fick
            for i = 1:height(t)
                RHVT = Geometry3D.RotMat2Fick(t.RDotinEyeRM{i});
                LHVT = Geometry3D.RotMat2Fick(t.LDotinEyeRM{i});

                t.RNewX(i) = rad2deg(RHVT(1));
                t.RNewY(i) = rad2deg(RHVT(2));
                t.RNewZ(i) = rad2deg(RHVT(3));
                t.LNewX(i) = rad2deg(LHVT(1));
                t.LNewY(i) = rad2deg(LHVT(2));
                t.LNewZ(i) = rad2deg(LHVT(3));
            end

            % Get disparity!!
            t.HDisparity = t.LNewX - t.RNewX;
            t.VDisparity = t.LNewY - t.RNewY;

            eyepoints = t;
        end


        function [screenPoints eyes] = PointsEyesToScreen(eyes, eyePoints, LeyeScreen, ReyeScreen)

            % TODO: correct for fick/helm or whatever coordinate systems we
            % have in the eye posints
            
            % TODO: Review how eye position is added (probably also needs
            % to be a rotation

            screenPoints.LX = LeyeScreen.middleX + LeyeScreen.distance*tand(eyes.L.H + eyePoints.LNewX)*LeyeScreen.pixPerCmWidth;
            screenPoints.LY = LeyeScreen.middleY + LeyeScreen.distance*tand(eyes.L.V + eyePoints.LNewY)*LeyeScreen.pixPerCmHeight;
            screenPoints.RX = ReyeScreen.middleX + ReyeScreen.distance*tand(eyes.R.H + eyePoints.RNewX)*ReyeScreen.pixPerCmWidth;
            screenPoints.RY = ReyeScreen.middleY + ReyeScreen.distance*tand(eyes.R.V + eyePoints.RNewY)*ReyeScreen.pixPerCmHeight;
                
        end

        function [heyes, hfix, hpoints, hLRpoints, hdisparity, hLscreen] = PlotPoints(points, eyePoints, screenPoints, eyes, fixationDistance)

            t = points;

            % View the dots, eyes, and fixation dot
            figure
            subplot(2,4,[1 2 5 6])

            hpoints=plot3(t.X,t.Z, t.Y,'o','Color','k'); hold on
            hfix=plot3(0,fixationDistance,0,'o','Color','r','LineWidth',2, 'markersize',20); % fixation spot
            heyes=plot3([eyes.L.X 0 eyes.R.X], [0 fixationDistance 0], [0 0 0],'-o','Color','b', 'markersize',20); % left eye fixation spot and right eye
            xlim(1.2*[min(t.X) max(t.X)])
            zlim(1.2*[min(t.Y) max(t.Y)])
            ylim([-5 min(100,max(fixationDistance, max(t.Z)))*1.2])
            xlabel('X (cm)')
            zlabel('Y (cm)')
            ylabel('Z (cm)')
            view(-30,30)

            hLRpoints = [];
            t = eyePoints;
            subplot(2,4,3)
            hLRpoints.R = plot(t.RNewX, t.RNewY,'o');
            hold
            hLRpoints.L = plot(t.LNewX, t.LNewY,'o');
            hLRpoints.FP = plot(0,0,'ro','linewidth',3);
            grid
            xlim([min([t.RNewX t.LNewX],[],"all") max([t.RNewX t.LNewX],[],"all")])
            ylim([min([t.RNewY t.LNewY],[],"all") max([t.RNewY t.LNewY],[],"all")])
            xlabel('X (deg)')
            ylabel('Y (deg)')
            legend({'Right Eye','Left Eye'})

            subplot(2,4,4)
            hdisparity = quiver(t.RNewX, t.RNewY, t.LNewX-t.RNewX, t.LNewY-t.RNewY);
            % for i = 1:height(t)
            %     line([t.RAnglex(i,1) t.LAnglex(i,1)],[t.RAngley(i,1) t.LAngley(i,1)])
            % end
            title('Disparity')
            set(gca,'xlim',[-20 20],'ylim',[-20 20])

            % subplot(2,4,7)
            % for i = 1:height(t)
            %     line([t.RNewX(i,1) t.LNewX(i,1)],[t.RNewY(i,1) t.LNewY(i,1)])
            % end
            % title('Disparity With Torsion')

            subplot(2,4,7)
            hLscreen.LPoints = plot(screenPoints.LX, screenPoints.LY, '+');
            hold
            hLscreen.LFP = plot(screenPoints.LX(end), screenPoints.LY(end), 'ro');
            set(gca,'xlim',[0 1920],'ylim',[0 1080])

            subplot(2,4,8)
            hLscreen.RPoints = plot(screenPoints.RX, screenPoints.RY, '+');
            hold
            hLscreen.RFP = plot(screenPoints.RX(end), screenPoints.RY(end), 'ro');
            set(gca,'xlim',[0 1920],'ylim',[0 1080])
            
        end

        function demoDisparityUpdate(sliderValues,heyes, hfix, hpoints,hLRpoints, hdisparity, hscreen, points)


            sizeCm = [30*16/9 30];
            resPix = [1920 1080];
            distance = 57;
            slant = 0;
            LeyeScreen = Geometry3D.MakeScreen(sizeCm, resPix, distance, slant);
            ReyeScreen = Geometry3D.MakeScreen(sizeCm, resPix, distance, slant);

            fixationSpot = [];
            fixationSpot.X = 0;
            fixationSpot.Y = 0;  
            fixationSpot.Z = sliderValues.fixationDistance;


            % points = GeneratePoints(sliderValues.nDots, sliderValues.sizeStimCm, sliderValues.stimDistance);
            eyes = Geometry3D.MakeEyes(sliderValues.ipdmm/10, fixationSpot, sliderValues.torsionVersion, sliderValues.torsionVergence);
            points.Z = ones(size(points.Z))*sliderValues.stimDistance;


            points{end,:} = [fixationSpot.X fixationSpot.Y fixationSpot.Z];

            [t] = Geometry3D.Points3DToEyes(points, eyes);
            screenPoints = Geometry3D.PointsEyesToScreen(eyes, t, LeyeScreen, ReyeScreen);

            set(hfix, 'ydata', sliderValues.fixationDistance);
            set(heyes, 'xdata', [eyes.L.X 0 eyes.R.X], 'ydata', [0 sliderValues.fixationDistance 0], 'zdata', [0 0 0]);
            set(hpoints, 'xdata', points.X, 'ydata',points.Z, 'zdata', points.Y);

            set(hLRpoints.L, 'xdata', t.LNewX, 'ydata', t.LNewY);
            set(hLRpoints.R, 'xdata', t.RNewX, 'ydata', t.RNewY);
            


            set(hdisparity, 'xdata',(t.LNewX+t.RNewX)/2, 'ydata', (t.LNewY+t.RNewY)/2, 'udata', t.LNewX-t.RNewX, 'vdata', t.LNewY-t.RNewY);

            set(hscreen.LPoints, 'xdata', screenPoints.LX, 'ydata',screenPoints.LY);
            set(hscreen.RPoints, 'xdata', screenPoints.RX, 'ydata',screenPoints.RY);            
            set(hscreen.LFP, 'xdata', screenPoints.LX(end), 'ydata', screenPoints.LY(end));            
            set(hscreen.RFP, 'xdata', screenPoints.RX(end), 'ydata', screenPoints.RY(end));  

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