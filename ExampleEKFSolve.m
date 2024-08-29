%%
clear all, close all


stim = "Ground plane"; % "Ground plane" or "Sphere at 1m"
w = [0 -0.1 0]'; %rad/s
v = [1 0 0]'; %m/s
eyeElevation = -30; %deg
h =  1.5; % m
N = 250;

noiseLevelStd = 0.0; 
Niter = 50;

% Calculate the true motion field
[motionField, visualDirections, motionFieldLinear, motionFieldRotational, Jv, Jw] = Geometry3D.CalculateMotionField(N, w, v, h, eyeElevation, stim);

% stack everything for the two dimensions of motion
visualDirections = [visualDirections;visualDirections];
motionField = [motionField(:,1);motionField(:,2)];
Jvv = [ squeeze(Jv(1,:,:))';squeeze(Jv(2,:,:))'];
Jww = [ squeeze(Jw(1,:,:))';squeeze(Jw(2,:,:))'];


% add noise to the motion field
motionFieldMeasurement = motionField + randn(size(motionField))*deg2rad(noiseLevelStd);


% true values of the estimates
trueDepth = max( -1/h * (visualDirections(1:end/2,1)*sind(eyeElevation) + visualDirections(1:end/2,3)*cosd(eyeElevation)) , 0);
trueState =[w;v;trueDepth];


% initial estimates
wEst = [0 0 0]';
vEst = [0 0 0]';
eyeelEst = 0;
hest = h; % known
% Dest = zeros(height(visualDirections)*2,1);
DepthEst = max( -1/hest * (visualDirections(:,1)*sind(-20) + visualDirections(:,3)*cosd(-20)) , 0);
% DepthEst = rand(size(x));
DepthEst = DepthEst(1:end/2);

xState = [wEst; vEst; DepthEst];
xStateError = zeros(numel(xState),Niter);

F = eye(height(xState));
Q = 0.1*diag([1 1 1 1 1 1 1*ones(1,height(DepthEst))]);% eye(height(xState))*.1; % state noise
R = eye(height(motionFieldMeasurement))*0.01; % measurement noise
% P = diag([1 1 1 1 1 1 10*ones(1,height(DepthEst))]); % State covariance estimate

% P = 0.1*[diag([1 1 1 1 1 1])  zeros(6, height(DepthEst)); zeros(height(DepthEst), 6) eye(height(DepthEst))]; % State covariance estimate

P = eye(height(xStateError));

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
    hEstimate = Jww*wEst + D*Jvv*vEst;
    yInnovation = motionFieldMeasurement - hEstimate;


    % Jacobian
    hw = Jww;
    hv = D*Jvv;
    hDd = Jvv*vEst;
    hD = [diag(0.5*(hDd(1:end/2) + hDd(end/2+1:end))); diag(0.5*(hDd(1:end/2) + hDd(end/2+1:end)))];


    H = [hw hv hD];

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

figure
subplot(3,2,1)
plot(xStateError(1:6,:)','-o')
hold
plot(mean(xStateError(1:6,:)),'r','linewidth',2)
ylabel('Error in state estimates')
xlabel('iteration')
legend({'wx' 'wy' 'wz' 'vx' 'vy' 'vz' 'Average'})
subplot(3,2,2)
plot(xStateError(7:end,:)','color',[0.6 0.6 0.6]);
hold
plot(mean(xStateError(7:end,:)),'r','linewidth',2)
% set(gca,'ylim',[-1 1])
ylabel('Error in depth estimates')
xlabel('iteration')


% plot motion fields
subplot(3,2,3 , 'nextplot' , 'add')
title('motion field (normalized)')
colors = get(gca,'colororder');
hs = [];

az = visualDirections(:,2).*acos(visualDirections(:,1));
el = visualDirections(:,3).*acos(visualDirections(:,1));
azdeg = rad2deg(az);
eldeg = rad2deg(el);
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
quiver( azdeg(1:end/2),  eldeg(1:end/2), zeros(size(trueDepth)), trueDepth,'linewidth',1.5 ) ;
quiver( azdeg(1:end/2),  eldeg(1:end/2), zeros(size(trueDepth)), DepthEst,'color','r');
legend({'True' 'Estimate'})


subplot(3,2,6 , 'nextplot' , 'add')
title('depth error (normalized)')
quiver( azdeg(1:end/2),  eldeg(1:end/2), zeros(size(trueDepth)), trueDepth- DepthEst,'color','r');