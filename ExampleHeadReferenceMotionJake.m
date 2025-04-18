%%
clear all, close all

%% Set parameters 

% x is positive forward
% y is positive right
% z is positive up

% heading velocity in head reference 
headingSpeed    = 10; % m/s
headingAzimuth  = 0; % deg postiive up
headingVelocity = headingSpeed*[0 0 -1 ]'; 

% eye position
eyePositionHeight   = 1.5; % m positive up
eyeAzimuth          = 30; % deg positive right
eyeElevation        = -90; % deg positive up

% eye movement gain
gain = 1;
% gain = 0.5;

plotLimit = 90;
numberOfPoints = 2000;


%% Calculate eye position and eye velocity
% get the rotation matrix of the eye (doing a listings rotation so there is
% no false torsion) 
direction = atan2(eyeElevation, -eyeAzimuth);
eccentricity = deg2rad(sqrt(eyeElevation.^2 + eyeAzimuth.^2));
eyePositionRotMat = Geometry3D.List2Mat([direction eccentricity 0]);
% get the eye location in the world
eyeLocation = [0, 0, eyePositionHeight];

% calculate the eye linear velocity by multiplying the velocity vector by
% the rotation matrix of the eye position 
eyeLinearVelocity = eyePositionRotMat'*headingVelocity;

% calculate the angular eye velocity assuming the eye is pursuing the point
% in the motion field caused by the heading velocity (times a gain factor)
foveaDirection = [1 0 0];
fovealVelocity = Geometry3D.CalculateMotionField(foveaDirection, [0 0 0]', eyeLinearVelocity, eyeLocation, eyePositionRotMat);
[~, JacAngular] = Geometry3D.CalculateMotionJacobianFields(foveaDirection);

eyeAngularVelocity = -gain*JacAngular'*fovealVelocity'; 

% todo: think about torsion. This formula and the way Jw was calculated
% implies pursuit following listing's law with axis in the ZY plane  


%% Calculate motion fields
N = numberOfPoints;
visualDirections = Geometry3D.SampleVisualDirections(N,'Random');

% 1 - head reference motion field (no eye rotation)
% 2 - eye reference motion field without eye velocity
% 3 - eye reference motion field with eye velocity
motionFieldLinearHeadReference = Geometry3D.CalculateMotionField(visualDirections, [0 0 0]', headingVelocity, eyeLocation); 
motionFieldLinearEyeReference = Geometry3D.CalculateMotionField(visualDirections, [0 0 0]', eyeLinearVelocity, eyeLocation, eyePositionRotMat);
motionFieldTotalEyeRefefence = Geometry3D.CalculateMotionField(visualDirections, eyeAngularVelocity, eyeLinearVelocity, eyeLocation, eyePositionRotMat);




%% Plot the motion fields

figure('Color','w', 'name','generation')
headingVelocity = headingVelocity ./ vecnorm(headingVelocity);
eyeLinearVelocity = eyeLinearVelocity ./ vecnorm(eyeLinearVelocity);

vazHead = cos(atan2(headingVelocity(3),-headingVelocity(2))).*acosd(headingVelocity(1));
velHead = 0;


vazEye = cos(atan2(eyeLinearVelocity(3),-eyeLinearVelocity(2))).*acosd(eyeLinearVelocity(1));
velEye = sin(atan2(eyeLinearVelocity(3),-eyeLinearVelocity(2))).*acosd(eyeLinearVelocity(1));

eyeAz = cos(atan2(eyePositionRotMat(3,1),-eyePositionRotMat(2,1))).*acosd(eyePositionRotMat(1,1));
eyeEl = sin(atan2(eyePositionRotMat(3,1),-eyePositionRotMat(2,1))).*acosd(eyePositionRotMat(1,1));

tiledlayout(1,3);

nexttile
hq1 = plotMotionField(visualDirections, motionFieldLinearHeadReference);

title({'Head reference' sprintf('(heading %0.1f m/s az. %0.1f deg)',headingSpeed, headingAzimuth)})
h1 = plot(vazHead,velHead,'go','linewidth',2);
h2 = plot(eyeAz,eyeEl,'ro','linewidth',2);
legend([h1,h2],{'heading' 'gaze'},'box','off','fontsize',14)
set(gca,'xlim',[-plotLimit plotLimit], 'ylim',[-plotLimit plotLimit])

nexttile
hq2 = plotMotionField(visualDirections, motionFieldLinearEyeReference);
title({'Eye reference' '(eye not moving)'})
h1 = plot(vazEye,velEye,'go','linewidth',2);
h2 = plot(0,0,'ro','linewidth',2);
legend([h1,h2],{'heading' 'gaze'},'box','off','fontsize',14)
set(gca,'xlim',[-plotLimit plotLimit], 'ylim',[-plotLimit plotLimit])

nexttile
hq3 = plotMotionField(visualDirections, motionFieldTotalEyeRefefence);
title({'Eye reference' sprintf('(eye moving with gain %0.1f)',gain)})
h1 = plot(vazEye,velEye,'go','linewidth',2);
h2 = plot(0,0,'ro','linewidth',2);
legend([h1,h2],{'heading' 'gaze'},'box','off','fontsize',14)
set(gca,'xlim',[-plotLimit plotLimit], 'ylim',[-plotLimit plotLimit])

function hq = plotMotionField(visualDirections, motionField)

set(gca,'nextplot','add');

%line([-90 90],[0 0],'linestyle','-.','color',0.5*[1 1 1])
%line([0 0],[-90 90],'linestyle','-.','color',0.5*[1 1 1])
t =0:1:360;x=cosd(t);y=sind(t);
plot(90*x,90*y,'linestyle','-.','color',0.5*[1 1 1])
plot(60*x,60*y,'linestyle','-.','color',0.5*[1 1 1])
plot(30*x,30*y,'linestyle','-.','color',0.5*[1 1 1])
axis square
%set(gca,'visible','off')
text(-30*cosd(45), 30*sind(45), '30','HorizontalAlignment','right','VerticalAlignment','bottom','color',0.5*[1 1 1])
text(-60*cosd(45), 60*sind(45), '60','HorizontalAlignment','right','VerticalAlignment','bottom','color',0.5*[1 1 1])
text(-90*cosd(45), 90*sind(45), '90','HorizontalAlignment','right','VerticalAlignment','bottom','color',0.5*[1 1 1])
set(gca,'xlim',[-90 90])
set(gca,'ylim',[-90 90])
set(gca,'xtick',[],'ytick',[],'xcolor',0.5*[1 1 1],'ycolor',0.5*[1 1 1]);
set(gca,'YAxisLocation','origin','XAxisLocation','origin')


% normalize fields for plotting (maybe)
motionField = motionField ./ vecnorm(motionField, 2, 2);
autoscale = 'on';
motionAz = cos(atan2(visualDirections(:,3),-visualDirections(:,2))).*acosd(visualDirections(:,1));
motionEl = sin(atan2(visualDirections(:,3),-visualDirections(:,2))).*acosd(visualDirections(:,1));

hq = quiver(motionAz,motionEl, motionField(:,1), -motionField(:,2),'color','k','linewidth',2, 'autoscale',autoscale);
end
