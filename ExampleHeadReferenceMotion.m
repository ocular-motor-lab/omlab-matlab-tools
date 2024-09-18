%%
clear all, close all
% heading velocity in head reference 
headingSpeed    = 10; % m/s
headingAzimuth  = 10; % deg postiive up
v = headingSpeed*[cosd(headingAzimuth) -sind(headingAzimuth) 0]'; 

% eye position
eyePositionHeight   = 1.5; % m positive up
eyeAzimuth          = 6.1; % deg positive right
eyeElevation        = -7.12; % deg positive up

% eye movement gain
gain = 1;
% gain = 0.5;

plotLimit = 40;
numberOfPoints = 5000;

% get the rotation matrix of the eye (doing a listings rotation so there is no false torsion)
direction = atan2(deg2rad(eyeElevation), -deg2rad(eyeAzimuth));
eccentricity = sqrt(deg2rad(eyeElevation).^2 + deg2rad(eyeAzimuth).^2);
Reye = Geometry3D.List2Mat([direction eccentricity 0]);

% calculate the linear velocity in eye coordinates
veye = Reye'*v;

% calculate the linear velocity at the gaze direction
linearV = Geometry3D.CalculateMotionField([1 0 0], [0 0 0]', veye, [0,0, eyePositionHeight], Reye);

% calculate the rotational eye velocity needed to cancel that motion (times
% a gain factor)
[Jv, Jw] = Geometry3D.CalculateMotionJacobianFields([1 0 0]);
w = -gain*Jw'*linearV'; % todo: think about torsion. This formula and the way Jw was calculated implies pursuit following listing's law with axis in the ZY plane 


% calculate motion fields
N = numberOfPoints;
visualDirections = Geometry3D.SampleVisualDirections(N);

% 1 - head reference motion field (no eye rotation)
% 2 - eye reference motion field without eye velocity
% 3 - eye reference motion field with eye velocity
motionField1 = Geometry3D.CalculateMotionField(visualDirections, [0 0 0]', v, [0,0, eyePositionHeight]); 
motionField2 = Geometry3D.CalculateMotionField(visualDirections, [0 0 0]', veye, [0,0, eyePositionHeight], Reye);
motionField3 = Geometry3D.CalculateMotionField(visualDirections, w, veye, [0,0, eyePositionHeight], Reye);

motionField1 = rad2deg(motionField1);
motionField2 = rad2deg(motionField2);
motionField3 = rad2deg(motionField3);

motionField1 = motionField1 ./ vecnorm(motionField1, 2, 2);
motionField2 = motionField2 ./ vecnorm(motionField2, 2, 2);
motionField3 = motionField3 ./ vecnorm(motionField3, 2, 2);

% plot the fields

figure('Color','w')
v = v ./ vecnorm(v);
veye = veye ./ vecnorm(veye);

vazHead = cos(atan2(v(3),-v(2))).*acosd(v(1));
velHead = 0;
vazEye = cos(atan2(veye(3),-veye(2))).*acosd(veye(1));
velEye = sin(atan2(veye(3),-veye(2))).*acosd(veye(1));

eyeAz = cos(atan2(Reye(3,1),-Reye(2,1))).*acosd(Reye(1,1));
eyeEl = sin(atan2(Reye(3,1),-Reye(2,1))).*acosd(Reye(1,1));

motionAz = cos(atan2(visualDirections(:,3),-visualDirections(:,2))).*acosd(visualDirections(:,1));
motionEl = sin(atan2(visualDirections(:,3),-visualDirections(:,2))).*acosd(visualDirections(:,1));

autoscale = 'on';

tiledlayout(1,3);
nexttile, set(gca,'nextplot','add');
AddAxes();
hq1 = quiver(motionAz,motionEl, motionField1(:,1), -motionField1(:,2),'color','k','linewidth',2, 'autoscale',autoscale);
title({'Head reference' sprintf('(heading %0.1f m/s az. %0.1f deg)',headingSpeed, headingAzimuth)})

h1 = plot(vazHead,velHead,'go','linewidth',2);
h2 = plot(eyeAz,eyeEl,'ro','linewidth',2);
legend([h1,h2],{'heading' 'gaze'},'box','off','fontsize',14)

set(gca,'xlim',[-plotLimit plotLimit], 'ylim',[-plotLimit plotLimit])

% add tv
t =0:1:360;
tvdeg = 60;
tvaspectratio = 16/9;
tvdegx = tvdeg/2;
tvdegy = tvdegx/tvaspectratio;
tvdiag = sqrt(tvdegx^2 + tvdegy.^2);

tvx = min(max(cosd(t)*tvdiag, -tvdegx),tvdegx);
tvy = min(max(sind(t)*tvdiag, -tvdegy),tvdegy);
tvz = ones(size(tvx)) * tvdegx/tand(tvdegx);
r = sqrt(tvx.^2 +tvy.^2 +tvz.^2);
tvx = tvx./r;
tvy = tvy./r;
tvz = tvz./r;

tvAz = cos(atan2(tvy,-tvx)).*acosd(tvz);
tvEl = sin(atan2(tvy,-tvx)).*acosd(tvz);

% figure
plot(tvAz, tvEl);
% a=1;


nexttile, set(gca,'nextplot','add');
AddAxes();
quiver(motionAz,motionEl, motionField2(:,1), -motionField2(:,2),'color','k','linewidth',2, 'autoscale',autoscale)
title({'Eye reference' '(eye not moving)'})
h1 = plot(vazEye,velEye,'go','linewidth',2);
h2 = plot(0,0,'ro','linewidth',2);
legend([h1,h2],{'heading' 'gaze'},'box','off','fontsize',14)
set(gca,'xlim',[-plotLimit plotLimit], 'ylim',[-plotLimit plotLimit])

nexttile, set(gca,'nextplot','add');
AddAxes();
quiver(motionAz,motionEl, motionField3(:,1), -motionField3(:,2),'color','k','linewidth',2, 'autoscale',autoscale)
title({'Eye reference' sprintf('(eye moving with gain %0.1f)',gain)})
h1 = plot(vazEye,velEye,'go','linewidth',2);
h2 = plot(0,0,'ro','linewidth',2);
legend([h1,h2],{'heading' 'gaze'},'box','off','fontsize',14)
set(gca,'xlim',[-plotLimit plotLimit], 'ylim',[-plotLimit plotLimit])


function AddAxes()
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

end