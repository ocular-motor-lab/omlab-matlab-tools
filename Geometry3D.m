classdef Geometry3D
    %Geometry3D Summary of this class goes here
    %   Detailed explanation goes here

    properties
    end

    methods(Static)

        function demoDisparity
            %%
            nDots = 90;
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

            eyes = MakeEyes(idpCm, fixationDistance, torsionVersion, torsionVergence);
            [eyePoints] = Points3DToEyes(points, eyes);
            [heyes, hfix, hpoints, hdisparity] = PlotPoints(points, eyePoints, eyes, fixationDistance);

            s = SliderUI('DisparityCalculator',@(sliderValues) (demoDisparityUpdate(sliderValues,heyes, hfix, hpoints,hdisparity, points)), 0.2);
            s.AddControl('fixationDistance',30, [10 200])
            s.AddControl('stimDistance',    40, [10 200])
            s.AddControl('ipdmm',           60, [10 100])
            s.AddControl('torsionVersion',  0,  [0 30])
            s.AddControl('torsionVergence', 0,  [0 10])
            s.AddControl('nDots',           90, [0 100])
            s.AddControl('sizeStimCm',      10, [0 100])
            s.Open();

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