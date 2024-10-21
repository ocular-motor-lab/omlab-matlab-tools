%%
% clear all, close all


stim = "Ground plane"; % "Ground plane" or "Sphere at 1m"
w = [0 0 0]'; %rad/s
v = [2 0 0]'; %m/s
eyeElevation = -20; %deg
h =  1.5; % m
N = 500;

noiseLevelStd = 0; 

[motionField, visualDirections, motionFieldLinear, motionFieldRotational, Jv, Jw] = Geometry3D.CalculateMotionField(N, w, v, h, eyeElevation, stim);

% add noise to the motion field
motionField = motionField + randn(size(motionField))*deg2rad(noiseLevelStd);



% stack everything
x = [visualDirections(:,1);visualDirections(:,1)];
y = [visualDirections(:,2);visualDirections(:,2)];
z = [visualDirections(:,3);visualDirections(:,3)];
motionField = [motionField(:,1);motionField(:,2)];
Jvv = [ squeeze(Jv(1,:,:))';squeeze(Jv(2,:,:))'];
Jww = [ squeeze(Jw(1,:,:))';squeeze(Jw(2,:,:))'];
J = [Jww,zeros(size(Jww));zeros(size(Jvv)),Jvv];
% some elements may be nan, not sure why ...
badSamples = isnan(motionField);
motionField(badSamples) = [];
x(badSamples) = [];
y(badSamples) = [];
z(badSamples) = [];
J(badSamples) = []; % REVIEW THIS


% true values of the estimates
truestate =[w;v];
trueD = max( -1/h * (x(1:end/2)*sind(eyeElevation) + z(1:end/2)*cosd(eyeElevation)) , 0);

% variables to collect the evolution of the errors
serror = [];
Derror = [];

% initial estimates
stateEst = [0 0 0 0 0 0]';
eyeelEst = 0;
hest = h; % known
% Dest = zeros(height(visualDirections)*2,1);
DepthEst = max( -1/hest * (x*sind(-20) + z*cosd(-20)) , 0);
% DepthEst = rand(size(x));
DepthEst = DepthEst(1:end/2);

% gradient descend
for i=1:500

    D = diag([DepthEst;DepthEst]);

    % sensory prediction error
    motionFieldPredictionError = (Jww*stateEst(1:3) + D*Jvv*stateEst(4:6) ) - motionField;

    % gradeints of loss of over the estimates
    % loss is mean squared motion field prediction error

    % assume we know that the pairs of motion estimates are matched to the
    % same visual direction and therefore the same depth. So we have half
    % the D estimates as motion predictions
    gradD = -2*motionFieldPredictionError.*Jvv*stateEst(4:6);
    gradD = (gradD(1:end/2) + gradD(end/2+1:end))/2;
    gradw = -2*Jww'*motionFieldPredictionError;
    gradv = -2*Jvv'*(motionFieldPredictionError.*[DepthEst;DepthEst]);

    % learning rates
    lr = 0.002;

    DepthEst = DepthEst + 10*lr*gradD;
    stateEst = stateEst + lr*[gradw;gradv];

    
    serror(:,i) = truestate - stateEst;
    Derror(:,i) = trueD - DepthEst;
end


figure
subplot(3,2,1)
plot(serror')
hold
plot(mean(serror),'r','linewidth',2)
ylabel('Error in state estimates')
xlabel('iteration')
legend({'wx' 'wy' 'wz' 'vx' 'vy' 'vz' 'Average'})
subplot(3,2,2)
plot(Derror','color',[0.6 0.6 0.6]);
hold
plot(mean(Derror),'r','linewidth',2)
% set(gca,'ylim',[-1 1])
ylabel('Error in depth estimates')
xlabel('iteration')


% plot motion fields
subplot(3,2,3 , 'nextplot' , 'add')
title('motion field (normalized)')
colors = get(gca,'colororder');
hs = [];

az = y.*acos(x);
el = z.*acos(x);
azdeg = rad2deg(az);
eldeg = rad2deg(el);
motionFieldRotational = rad2deg(motionFieldRotational);
motionFieldLinear = rad2deg(motionFieldLinear);
motionField = rad2deg(motionField);

quiver( azdeg(1:end/2),  eldeg(1:end/2), motionField(1:end/2), motionField(end/2+1:end),'linewidth',1.5 ) ;
    motionFieldPred = (Jww*stateEst(1:3) + D*Jvv*stateEst(4:6)) ;
quiver( azdeg(1:end/2),  eldeg(1:end/2), motionFieldPred(1:end/2), motionFieldPred(end/2+1:end),'color','r');
legend({'True' 'Prediction'})

subplot(3,2,4 , 'nextplot' , 'add')
title('motion field error (normalized)')
quiver( azdeg(1:end/2),  eldeg(1:end/2), motionFieldPred(1:end/2)- motionField(1:end/2), motionFieldPred(end/2+1:end)-motionField(end/2+1:end),'color','r');


subplot(3,2,5 , 'nextplot' , 'add')
title('depth (normalized up means closer)')
quiver( azdeg(1:end/2),  eldeg(1:end/2), zeros(size(trueD)), trueD,'linewidth',1.5 ) ;
quiver( azdeg(1:end/2),  eldeg(1:end/2), zeros(size(trueD)), DepthEst,'color','r');
legend({'True' 'Estimate'})


subplot(3,2,6 , 'nextplot' , 'add')
title('depth error (normalized)')
quiver( azdeg(1:end/2),  eldeg(1:end/2), zeros(size(trueD)), trueD- DepthEst,'color','r');