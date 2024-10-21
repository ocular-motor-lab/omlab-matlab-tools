%%
clear all, close all

w = [0 0 0]';
v = [1 0.5 0]';
eyeElevation = 0; %deg
h =  1;
N = 250;

stim = "Ground plane"; % "Ground plane" or "Sphere at 1m"

[motionField, visualDirections, motionFieldLinear, motionFieldRotational, Jv, Jw] = Geometry3D.CalculateMotionField(N, w, v, h, eyeElevation, stim);

motionField = motionField + randn(size(motionField))*2;

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





%% noise level simulations
clear all, close all

w = [0 0 0]';
% v = [1 0.5 0]';
v = [0 0 0]';
eyeElevation = 0; %deg
h =  1;
N = 250;

stim = "Ground plane"; % "Ground plane" or "Sphere at 1m"

[motionField, visualDirections, motionFieldLinear, motionFieldRotational, Jv, Jw] = Geometry3D.CalculateMotionField(N, w, v, h, eyeElevation, stim);
motionFieldo = motionField;

noises = 0:1:20;
reps = 300;
errors = zeros(6,length(noises), reps);

for i=1:length(noises)
    i
    for j=1:reps
        motionField = motionFieldo + randn(size(motionFieldo))*noises(i);


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
        west = state(1:3);
        vest = state(4:6);


        truestate =[w;v];
        error = truestate - state;

        errors(:,i,j) = error.*error;

    end
end

%%
figure
plot(noises, mean(errors,3),'linewidth',2)
legend({'wx','wy','wz','vx','vy','vz'},'fontsize',15)
xlabel('noise level')
ylabel('mean square error')


%%
clear all, close all

w = [10 0 0]';
v = [3 0.5 0]';
eyeElevation = 20; %deg
h =  1;
N = 250;

stim = "Ground plane"; % "Ground plane" or "Sphere at 1m"

[motionField, visualDirections, motionFieldLinear, motionFieldRotational, Jv, Jw] = Geometry3D.CalculateMotionField(N, w, v, h, eyeElevation, stim);

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


motionField = motionField + randn(size(motionField))*0;
serror = [];
Derror = [];
stateEst = [0 0 0 0 0 0]';
% vest = [0 0 0];
eyeelEst = 0;
hest = h; % known
truestate =[w;v];
trueD = max( -1/hest * (x*sind(eyeElevation) + z*cosd(eyeElevation)) , 0);
% Dest = zeros(height(visualDirections)*2,1);
Dest = max( -1/hest * (x*sind(0) + z*cosd(0)) , 0);
Dest = Dest(1:end/2);
trueD = trueD(1:end/2);

% gradient descend
for i=1:500

    D = diag([Dest;Dest]);

    % sensory prediction error
    motionFieldPredictionError = motionField - [eye(size(D)) D]*stateEst;

    % gradeints over the estimates
    gradD = motionFieldPredictionError(1:end/2).*Jvv(1:end/2,:)*stateEst(4:6)+2*motionFieldPredictionError(1:end/2).*Jvv(1:end/2,:)*stateEst(4:6);
    gradw = 2*Jww'*motionFieldPredictionError;
    gradv = 2*Jvv'*(motionFieldPredictionError.*[Dest;Dest]);

    % learning rates
    lr = 0.002;

    Dest = Dest + lr*gradD;
    stateEst = stateEst + lr*[gradw;gradv];

    
    serror(:,i) = truestate - stateEst;
    Derror(:,i) = trueD - Dest;
end

figure
subplot(2,1,1)
plot(serror','color',[0.6 0.6 0.6])
hold
plot(mean(serror),'r','linewidth',2)
ylabel('Error in state estimates')
subplot(2,1,2)
plot(Derror','color',[0.6 0.6 0.6])
hold
plot(mean(Derror),'r','linewidth',2)
ylabel('Error in depth estimates')
xlabel('iteration')
