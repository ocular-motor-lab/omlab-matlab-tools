%% Vergence

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

