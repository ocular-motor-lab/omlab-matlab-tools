%%
% clear all, close all
% heading velocity in head reference 
headingSpeed = 1; % m/s
headingAzimuth = 10; % deg postiive up
v = headingSpeed*[cosd(headingAzimuth) -sind(headingAzimuth) 0]'; 

% eye position
eyePositionHeight = 1.5; % m positive up
eyeAzimuth = 30; % deg positive right
eyeElevation = -40; % deg positive up

% eye movement gain
gain = 1;

% get the rotation matrix of the eye (TODO: think about torsion)
Reye = Geometry3D.Helm2RotMat([-deg2rad(eyeAzimuth),  -deg2rad(eyeElevation), 0]);

% calculate the linear velocity in eye coordinates
veye = Reye'*v;

% calculate the linear velocity at the gaze direction
linearV = Geometry3D.CalculateMotionField([1 0 0], [0 0 0]', veye, [0,0, eyePositionHeight], Reye);

% calculate the rotational eye velocity needed to cancel that motion (times
% a gain factor)
[Jv, Jw] = Geometry3D.CalculateMotionJacobianFields([1 0 0]);
w = -gain*Jw'*linearV';


% calculate motion fields
N = 250;
visualDirections = Geometry3D.SampleVisualDirections(N);

% 1 - head reference motion field (no eye rotation)
% 2 - eye reference motion field without eye velocity
% 3 - eye reference motion field with eye velocity

motionField1 = Geometry3D.CalculateMotionField(visualDirections, [0 0 0]', v, [0,0, eyePositionHeight]); 
motionField2 = Geometry3D.CalculateMotionField(visualDirections, [0 0 0]', veye, [0,0, eyePositionHeight], Reye);
motionField3 = Geometry3D.CalculateMotionField(visualDirections, w, veye, [0,0, eyePositionHeight], Reye);




figure

vaz = cos(atan2(v(3),-v(2))).*acosd(v(1));
vazeye = cos(atan2(veye(3),-veye(2))).*acosd(veye(1));
veleye = sin(atan2(veye(3),-veye(2))).*acosd(veye(1));

az = cos(atan2(visualDirections(:,3),-visualDirections(:,2))).*acosd(visualDirections(:,1));
el = sin(atan2(visualDirections(:,3),-visualDirections(:,2))).*acosd(visualDirections(:,1));


subplot(1,3,1);
quiver(az,el, motionField1(:,1), -motionField1(:,2))
line([-90 90],[0 0])
line([0 0],[-90 90])
title(sprintf('Head reference (heading %0.1f m/s %az=%0.1f deg',headingSpeed, headingAzimuth))
hold
h1 = plot(vaz,0,'ro');
h2 = plot(eyeAzimuth,eyeElevation,'ko');
legend([h1,h2],{'heading' 'gaze'})

subplot(1,3,2); 
quiver(az,el, motionField2(:,1), -motionField2(:,2))
line([-90 90],[0 0])
line([0 0],[-90 90])
hold
title('Eye reference (eye not moving )')
h1 = plot(vazeye,veleye,'ro');
h2 = plot(0,0,'ko');
legend([h1,h2],{'heading' 'gaze'})

subplot(1,3,3); 
quiver(az,el, motionField3(:,1), -motionField3(:,2))
line([-90 90],[0 0])
line([0 0],[-90 90])
hold
title(sprintf('Eye reference (eye moving with gain %0.1f)',gain))
h1 = plot(vazeye,veleye,'ro');
h2 = plot(0,0,'ko');
legend([h1,h2],{'heading' 'gaze'})