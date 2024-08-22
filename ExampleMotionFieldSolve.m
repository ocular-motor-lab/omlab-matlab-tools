%%
clear all, close all

w = [0 0 0]';
v = [1 0.5 0]';
eyeElevation = 0; %deg
h =  1;
N = 250;

stim = "Ground plane"; % "Ground plane" or "Sphere at 1m"

[motionField, visualDirections, motionFieldLinear, motionFieldRotational, Jv, Jw] = Geometry3D.CalculateMotionField(N, w, v, h, eyeElevation, stim);

motionField = motionField + randn(size(motionField))*0.5;

west = [0 0 0];
vest = [0 0 0];
eyeelEst = 0;
hest = h; % known

% stack everything
x = [visualDirections(:,1);visualDirections(:,1)];
y = [visualDirections(:,2);visualDirections(:,2)];
z = [visualDirections(:,3);visualDirections(:,3)];
motionField = [motionField(:,1);motionField(:,2)];
Jvv = [ squeeze(Jv(1,:,:))';squeeze(Jv(2,:,:))'];
Jww = [ squeeze(Jw(1,:,:))';squeeze(Jw(2,:,:))'];
J = [Jww,zeros(size(Jww));zeros(size(Jvv)),Jvv];

badSamples = isnan(motionField);
motionField(badSamples) = []; 
x(badSamples) = []; 
y(badSamples) = []; 
z(badSamples) = []; 
J(badSamples) = []; % REVIEW THIS

D = diag( max( -1/hest * (x*sind(eyeElevation) + z*cosd(eyeElevation)) , 0)  );

ID = [eye(size(D)) D];
IDJ = ID*J;

state = IDJ\motionField;
west = state(1:3)
vest = state(4:6)
