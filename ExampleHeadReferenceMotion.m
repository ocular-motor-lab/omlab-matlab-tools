%%
clear all, close all

%% Set parameters 

% heading velocity in head reference 
headingSpeed    = 10; % m/s
headingAzimuth  = 10; % deg postiive up
headingVelocity = headingSpeed*[cosd(headingAzimuth) -sind(headingAzimuth) 0]'; 

% eye position
eyePositionHeight   = 1.5; % m positive up
eyeAzimuth          = 16.1; % deg positive right
eyeElevation        = -17.12; % deg positive up

% eye movement gain
gain = 1;
%= gain = 0.5;

plotLimit = 40;
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
visualDirections = Geometry3D.SampleVisualDirections(N);

% 1 - head reference motion field (no eye rotation)
% 2 - eye reference motion field without eye velocity
% 3 - eye reference motion field with eye velocity
motionFieldLinearHeadReference = Geometry3D.CalculateMotionField(visualDirections, [0 0 0]', headingVelocity, eyeLocation); 
motionFieldLinearEyeReference = Geometry3D.CalculateMotionField(visualDirections, [0 0 0]', eyeLinearVelocity, eyeLocation, eyePositionRotMat);
motionFieldTotalEyeRefefence = Geometry3D.CalculateMotionField(visualDirections, eyeAngularVelocity, eyeLinearVelocity, eyeLocation, eyePositionRotMat);

%% Solve to estimate heading only given some efference copy of eye angular velocity

% get an estimate of eye velocity (this could be coming from a mix of
% efference copy but also motion estimation). For now let's test what
% happens changing the gain and adding some noise
estimatedEyeVelocityGain = 0.5; 
estimatedEyeVelocityNoise = 0;
estimatedEyeAngularVelocity = eyeAngularVelocity*estimatedEyeVelocityGain + randn(size(eyeAngularVelocity))*estimatedEyeVelocityNoise;

% estimate teh motion field due to that estimated eye velocity and
% substract it from the full retinal field to estimate the motion field due
% to only linear motion
[estimatedAngularMotionField, ~, ~, JacLinear, JacAngular, DepthField] = Geometry3D.CalculateMotionField(visualDirections, estimatedEyeAngularVelocity, [0 0 0]', eyeLocation, eyePositionRotMat);
estimatedLinearMotionField = motionFieldTotalEyeRefefence - estimatedAngularMotionField;

% stack the horizontal and vertical measurements and jacobians to do the
% regression
JacLinearStacked = [ squeeze(JacLinear(1,:,:))';squeeze(JacLinear(2,:,:))'];
estimatedLinearMotionFieldStacked = [estimatedLinearMotionField(:,1);estimatedLinearMotionField(:,2)];
DepthFieldStacked = [DepthField zeros(size(DepthField));zeros(size(DepthField)) DepthField];
% solve via regression
estimatedLinearVelocity = (DepthFieldStacked*JacLinearStacked) \ estimatedLinearMotionFieldStacked;

estimatedHeadingVelocity = eyePositionRotMat*estimatedLinearVelocity;
trueHeadingVelocity = headingVelocity;

resultsOnlyLinearEstimation = table(estimatedHeadingVelocity, trueHeadingVelocity,  'VariableNames', {'estimatedHeadingVelocity', 'trueHeadingVelocity'})


%% Solve to estimate heading and angular velocity, both from motion field, no efference copy
% stack the horizontal and vertical measurements and jacobians to do the
% regression motion = [I D]*[Jang 0; 0; Jline][w v]
JacLinearStacked = [ squeeze(JacLinear(1,:,:))';squeeze(JacLinear(2,:,:))'];
JacAngularStacked = [ squeeze(JacAngular(1,:,:))';squeeze(JacAngular(2,:,:))'];
Jac = [JacAngularStacked zeros(size(JacAngularStacked)); zeros(size(JacAngularStacked)) JacLinearStacked];
motionFieldTotalEyeRefefenceStacked = [motionFieldTotalEyeRefefence(:,1);motionFieldTotalEyeRefefence(:,2)];

% We need a factor that depends on depth (in diopters) for the linear field
% but a factor that is just zero or one depending on the presence or
% absence of stimuli for the angular field. In this case only the floor 
DepthFieldStacked1 = [DepthField zeros(size(DepthField));zeros(size(DepthField)) DepthField];
DepthFieldStacked = [diag(diag(DepthFieldStacked1)>0) DepthFieldStacked1];

% solve via regression
estimatedAngularLinearVelocity = (DepthFieldStacked*Jac) \ motionFieldTotalEyeRefefenceStacked;

estimatedHeadingVelocity = eyePositionRotMat*estimatedAngularLinearVelocity(4:6);
estimatedEyeVelocity = estimatedAngularLinearVelocity(1:3);
trueHeadingVelocity = headingVelocity;
trueEyeVelocity = eyeAngularVelocity;

resultsFullEstimation = table(estimatedHeadingVelocity, trueHeadingVelocity, estimatedEyeVelocity, trueEyeVelocity, 'VariableNames', {'estimatedHeadingVelocity', 'trueHeadingVelocity', 'estimatedEyeVelocity', 'trueEyeVelocity'})


%% Solve to estimate heading using efference copy as a prior and having uniform noise across the field

efferencePriorSigma = 10;
headingPriorSigma = 10;
efferencePriorBiasPercentGain = 0;
measurementNoiseSigma = 0.0001;

% stack the horizontal and vertical measurements and jacobians to do the
% regression motion = [I D]*[Jang 0; 0; Jline][w v]
JacLinearStacked = [ squeeze(JacLinear(1,:,:))';squeeze(JacLinear(2,:,:))'];
JacAngularStacked = [ squeeze(JacAngular(1,:,:))';squeeze(JacAngular(2,:,:))'];
Jac = [JacAngularStacked zeros(size(JacAngularStacked)); zeros(size(JacAngularStacked)) JacLinearStacked];
motionFieldTotalEyeRefefenceStacked = [motionFieldTotalEyeRefefence(:,1);motionFieldTotalEyeRefefence(:,2)];

% We need a factor that depends on depth (in diopters) for the linear field
% but a factor that is just zero or one depending on the presence or
% absence of stimuli for the angular field. In this case only the floor 
DepthFieldStacked1 = [DepthField zeros(size(DepthField));zeros(size(DepthField)) DepthField];
DepthFieldStacked = [diag(diag(DepthFieldStacked1)>0) DepthFieldStacked1];


% the prior for heading is zero, the prior for angular velocity comes from
% the efference copy
MuPrior = [eyeAngularVelocity; 0; 0; 0];
SigmaPrior = diag([ efferencePriorSigma*ones(1,3) headingPriorSigma*ones(1,3)]);
SigmaLikelihood = diag(measurementNoiseSigma*ones(1,height(motionFieldTotalEyeRefefenceStacked)));
% Q: should we rotate the prior from eye reference to head reference? Right
% now it's in eye reference
% Q: how worse would it get if we use a prior for depth instead of the true
% depth?

% solve via WLS Ridge regression

X = DepthFieldStacked*Jac; % design matrix
y = motionFieldTotalEyeRefefenceStacked + randn(size(motionFieldTotalEyeRefefenceStacked))*measurementNoiseSigma; % add noise to measurements

estimatedAngularLinearVelocity = (X'/SigmaLikelihood*X + eye(length(MuPrior))/SigmaPrior) \ (X'/SigmaLikelihood*y + SigmaPrior\MuPrior);
estimatedAngularLinearVelocityCoVar= (X'/SigmaLikelihood*X + eye(length(MuPrior))/SigmaPrior) \ X'/SigmaLikelihood*X  /(X'/SigmaLikelihood*X + eye(length(MuPrior))/SigmaPrior);

estimatedHeadingVelocity = eyePositionRotMat*estimatedAngularLinearVelocity(4:6);
estimatedEyeVelocity = estimatedAngularLinearVelocity(1:3);
trueHeadingVelocity = headingVelocity;
trueEyeVelocity = eyeAngularVelocity;

resultsFullEstimationWithNoiseAndPriorEfference = table(estimatedHeadingVelocity, trueHeadingVelocity, (trueHeadingVelocity-estimatedHeadingVelocity)./trueHeadingVelocity*100, ...
    estimatedEyeVelocity, trueEyeVelocity, (trueEyeVelocity-estimatedEyeVelocity)./trueEyeVelocity*100, ...
    'VariableNames', {'estimatedHeadingVelocity', 'trueHeadingVelocity', 'heading error %', 'estimatedEyeVelocity', 'trueEyeVelocity','angular error %'})


%% 
if (0)

    xState = MuPrior;

    Q = 0.1*diag([1 1 1 1 1 1 1*ones(1,height(DepthEst))]);% eye(height(xState))*.1; % state noise
    R = SigmaLikelihood;
    % P = diag([1 1 1 1 1 1 10*ones(1,height(DepthEst))]); % State covariance estimate

    % gradient descend
    for i=1:Niter

        xStateError(:,i) = trueState - xState;


        wEst = xState(1:3);
        vEst = xState(4:6);
        DepthEst = xState(7:end);

        % assume we know that the pairs of motion estimates are matched to the
        % same visual direction and therefore the same depth. So we have half
        % the D estimates as motion predictions
        D = diag([DepthEst;DepthEst]);


        % state update just assumes the same state.
        % TODO: add dependency of position on velocities
        xState = xState;

        P = F*P*F' + Q;

        % sensory prediction error or innovation
        hEstimate = X*xState;
        yInnovation = y - hEstimate;


        H = X;

        % covariance of the innovation
        S = H*P*H' + R;

        % S = S + 1*eye(size(S));

        % kalman gain (Gradient)
        K = P*H'/S;

        % State update
        xState = xState + K*yInnovation;

        % update of covariance estimate
        P = (eye(height(P)) - K*H)*P;


    end
end



%% Plot the motion fields

figure('Color','w', 'name','generation')
headingVelocity = headingVelocity ./ vecnorm(headingVelocity);
eyeLinearVelocity = eyeLinearVelocity ./ vecnorm(eyeLinearVelocity);
estimatedHeadingVelocity = estimatedHeadingVelocity ./ vecnorm(estimatedHeadingVelocity);
estimatedLinearVelocity = estimatedLinearVelocity ./ vecnorm(estimatedLinearVelocity);

vazHead = cos(atan2(headingVelocity(3),-headingVelocity(2))).*acosd(headingVelocity(1));
velHead = 0;


estVazHead = cos(atan2(estimatedHeadingVelocity(3),-estimatedHeadingVelocity(2))).*acosd(estimatedHeadingVelocity(1));
estVelHead = sin(atan2(estimatedHeadingVelocity(3),-estimatedHeadingVelocity(2))).*acosd(estimatedHeadingVelocity(1));
estVazEye = cos(atan2(estimatedLinearVelocity(3),-estimatedLinearVelocity(2))).*acosd(estimatedLinearVelocity(1));
estVelEye = sin(atan2(estimatedLinearVelocity(3),-estimatedLinearVelocity(2))).*acosd(estimatedLinearVelocity(1));


vazEye = cos(atan2(eyeLinearVelocity(3),-eyeLinearVelocity(2))).*acosd(eyeLinearVelocity(1));
velEye = sin(atan2(eyeLinearVelocity(3),-eyeLinearVelocity(2))).*acosd(eyeLinearVelocity(1));

eyeAz = cos(atan2(eyePositionRotMat(3,1),-eyePositionRotMat(2,1))).*acosd(eyePositionRotMat(1,1));
eyeEl = sin(atan2(eyePositionRotMat(3,1),-eyePositionRotMat(2,1))).*acosd(eyePositionRotMat(1,1));

tiledlayout(2,3);

nexttile
hq1 = plotMotionField(visualDirections, motionFieldLinearHeadReference);

title({'Head reference' sprintf('(heading %0.1f m/s az. %0.1f deg)',headingSpeed, headingAzimuth)})
h1 = plot(vazHead,velHead,'go','linewidth',2);
h2 = plot(eyeAz,eyeEl,'ro','linewidth',2);
h3 = plot(estVazHead,estVelHead,'bo','linewidth',2);
AddTVframe();
legend([h1,h2,h3],{'heading' 'gaze' 'estimated heading'},'box','off','fontsize',14)
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


nexttile
hq4 = plotMotionField(visualDirections, estimatedLinearMotionField);
title({'Estimated linear field in retina' sprintf('(eff copy gain %0.1f)',estimatedEyeVelocityGain)})
h1 = plot(estVazEye,estVelEye,'bo','linewidth',2);
h2 = plot(0,0,'ro','linewidth',2);
legend([h1,h2],{'estimated heading' 'gaze'},'box','off','fontsize',14)
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



function AddTVframe()
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
end