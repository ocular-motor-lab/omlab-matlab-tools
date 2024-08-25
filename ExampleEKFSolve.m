%%
clear all, close all


stim = "Ground plane"; % "Ground plane" or "Sphere at 1m"
w = [0 0 0]'; %rad/s
v = [1 0.5 0]'; %m/s
eyeElevation = -90; %deg
h =  1.5; % m
N = 200;

noiseLevelStd = 0.5; 

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


% true values of the estimates
truestate =[w;v];
trueD = max( -1/h * (x(1:end/2)*sind(eyeElevation) + z(1:end/2)*cosd(eyeElevation)) , 0);

% variables to collect the evolution of the errors
serror = [];
Derror = [];

% initial estimates
wEst = [0 0 0]';
vEst = [0 0 0]';
eyeelEst = 0;
hest = h; % known
% Dest = zeros(height(visualDirections)*2,1);
DepthEst = max( -1/hest * (x*sind(-20) + z*cosd(-20)) , 0);
% DepthEst = rand(size(x));
DepthEst = DepthEst(1:end/2);

xState = [wEst; vEst; DepthEst];

F = zeros(height(xState));
Q = eye(height(xState))*0.01; % state noise
R = eye(height(motionField))*0.1; % measurement noise
P = eye(height(xState))*1; % State covariance estimate

% gradient descend
for i=1:100

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
    hEstimate = Jww*wEst + D*Jvv*vEst;
    yInnovation = motionField - hEstimate;


    % Jacobian
    hw = Jww;
    hv = D*Jvv;
    hDd = Jvv*vEst;
    hD = [diag(0.5*(hDd(1:end/2) + hDd(end/2+1:end))); diag(0.5*(hDd(1:end/2) + hDd(end/2+1:end)))];


    H = [hw hv hD];

    % covariance of the innovation
    S = H*P*H' + R;

    % kalman gain (Gradient)
    K = P*H'/S;

    % State update
    xState = xState + K*yInnovation;

    % update of covariance estimate
    P = (eye(height(P)) - K*H)*P;

    
    serror(:,i) = truestate - xState(1:6);
    Derror(:,i) = trueD - xState(7:end);
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
    motionFieldPred = (Jww*xState(1:3) + D*Jvv*xState(4:6)) ;
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