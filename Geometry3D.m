classdef Geometry3D
    %Geometry3D Summary of this class goes here
    %   Detailed explanation goes here

    properties
    end

    methods(Static)

        function demoDisparity
            %%
            numDots = 90;
            sizeStimCm = 20;
            stimDistance = 10;


            % Make some dots
            % Make a table to start putting everything together
            points = table();
            points.X=rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
            points.Y=rand(numDots,1)*sizeStimCm - (sizeStimCm/2);
            points.Z=ones(numDots, 1)*stimDistance;

            idpCm = 6;
            fixationDistance = 400;
            torsionVersion = 30;
            torsionVergence = 0;

            eyes = Geometry3D.MakeEyes(idpCm, fixationDistance, torsionVersion, torsionVergence);
            [eyePoints] = Geometry3D.Points3DToEyes(points, eyes);
            [heyes, hfix, hpoints, hdisparity] = Geometry3D.PlotPoints(points, eyePoints, eyes, fixationDistance);

            s = InteractiveUI('DisparityCalculator',@(sliderValues) (demoDisparityUpdate(sliderValues,heyes, hfix, hpoints,hdisparity, points)), 0.2);
            s.AddControl('fixationDistance',30, [10 200])
            s.AddControl('stimDistance',    40, [10 200])
            s.AddControl('ipdmm',           60, [10 100])
            s.AddControl('torsionVersion',  0,  [0 30])
            s.AddControl('torsionVergence', 0,  [0 10])
            s.AddControl('nDots',           90, [0 100])
            s.AddControl('sizeStimCm',      10, [0 100])
            s.Open();

        end

        function eyes = MakeEyes(ipdCm, fixationDistance, torsionVersion, torsionVergence)

            % Position of eyes
            eyes.R.X =  ipdCm/2;
            eyes.R.Y = 0;
            eyes.R.Z = 0;
            eyes.L.X = -ipdCm/2;%
            eyes.L.Y = 0;
            eyes.L.Y = 0;

            % Orientation of eyes
            eyes.R.H = -atan2d(eyes.R.X, fixationDistance);
            eyes.R.V = 0;
            eyes.R.T = torsionVersion - torsionVergence;
            eyes.L.H = -atan2d(eyes.L.X, fixationDistance);
            eyes.L.V = 0;
            eyes.L.T = torsionVersion + torsionVergence;

        end


        function [eyepoints eyes] = Points3DToEyes(points, eyes )


            t =table();


            % Get the new xs in cm for the two eyes (transform the dots from head reference
            % to eye reference [nothing about eye angle here!])
            t.newxR = points.X - eyes.R.X; % both in cm
            t.newxL = points.X - eyes.L.X; % no need to do this for y bc the two eyes are only displaced horizontally on the head

            % Get the ANGLES of the dots relative to the (straight ahead) eye position in cm
            t.DotsAngleRx = atan2d(t.newxR, points.Z);
            t.DotsAngleLx = atan2d(t.newxL, points.Z);
            t.DotsAngleRy = atan2d(points.Y, (points.Z)./abs(cosd(t.DotsAngleRx)));
            t.DotsAngleLy = atan2d(points.Y, (points.Z)./abs(cosd(t.DotsAngleLx)));

            % Calculate the angles of the dots in the eye reference frame accounting
            % for the angles (vergence) of the eyes
            t.RAnglex = t.DotsAngleRx - eyes.R.H; % both in deg
            t.LAnglex = t.DotsAngleLx - eyes.L.H;
            t.RAngley = t.DotsAngleRy;
            t.LAngley = t.DotsAngleLy;

            % Get disparity!!
            t.HDisparity = t.LAnglex - t.RAnglex;
            t.VDisparity = t.LAngley - t.RAngley;

            % Get the rotation matrices describing the eye angles
            REyeRMx = Geometry3D.RotZ(deg2rad(eyes.R.H)); % all rows are the same so we can just use the first one
            REyeRMy = Geometry3D.RotY(deg2rad(eyes.R.V));
            REyeRMz = Geometry3D.RotX(deg2rad(eyes.R.T));
            REyeRM = REyeRMx*REyeRMy*REyeRMz; % Rot matrix defining the orientation of the eye relative to straight ahead direction.

            LEyeRMx = Geometry3D.RotZ(deg2rad(eyes.L.H)); % all rows are the same so we can just use the first one
            LEyeRMy = Geometry3D.RotY(deg2rad(eyes.L.V));
            LEyeRMz = Geometry3D.RotX(deg2rad(eyes.L.T));
            LEyeRM = LEyeRMx*LEyeRMy*LEyeRMz; % Rot matrix defining the orientation of the eye relative to straight ahead direction.

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

            eyepoints = t;
        end

        function [heyes, hfix, hpoints, hdisparity] = PlotPoints(points, eyePoints, eyes, fixationDistance)

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

            t = eyePoints;
            subplot(2,4,3)
            plot(t.RAnglex, t.RAngley,'o')
            hold
            plot(t.LAnglex, t.LAngley,'o')
            plot(0,0,'ro','linewidth',3);
            grid
            xlim([min([t.RAnglex t.LAnglex],[],"all") max([t.RAnglex t.LAnglex],[],"all")])
            ylim([min([t.RAngley t.LAngley],[],"all") max([t.RAngley t.LAngley],[],"all")])
            xlabel('X (deg)')
            ylabel('Y (deg)')
            legend({'Right Eye','Left Eye'})

            subplot(2,4,4)
            hdisparity = quiver(t.RAnglex, t.RAngley, t.LAnglex-t.RAnglex, t.LAngley-t.RAngley);
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
        end

        function demoDisparityUpdate(sliderValues,heyes, hfix, hpoints, hdisparity, points)


            % points = GeneratePoints(sliderValues.nDots, sliderValues.sizeStimCm, sliderValues.stimDistance);
            eyes = MakeEyes(sliderValues.ipdmm/10, sliderValues.fixationDistance, sliderValues.torsionVersion, sliderValues.torsionVergence);
            points.Z = ones(size(points.Z))*sliderValues.stimDistance;
            [t] = Points3DToEyes(points, eyes);

            set(hfix, 'ydata', sliderValues.fixationDistance);
            set(heyes, 'xdata', [eyes.L.X 0 eyes.R.X], 'ydata', [0 sliderValues.fixationDistance 0], 'zdata', [0 0 0]);
            set(hpoints, 'xdata', points.X, 'ydata',points.Z, 'zdata', points.Y);

            set(hdisparity, 'xdata',(t.LNewX+t.RNewX)/2, 'ydata', (t.LNewY+t.RNewY)/2, 'udata', t.LNewX-t.RNewX, 'vdata', t.LNewY-t.RNewY);

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