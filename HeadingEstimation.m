classdef HeadingEstimation


    methods(Static)
        function RunExampleTrialSimulation()

            headingSpeed    = 10; % m/s
            headingAzimuth  = 10; % deg postiive up

            eyePositionHeight   = 1.5; % m positive up
            eyeAzimuth          = 16.1; % deg positive right
            eyeElevation        = -17.12; % deg positive up

            eyeMovementGain = 1;

            numberOfRetinalPoints = 2000;
            screenDistance = 1;
            screenWidth = tand(30)*screenDistance*2;
            screenAspectRatio = 16/9;

            params.expParams = HeadingEstimation.SetUpExperiment(headingSpeed, headingAzimuth, eyePositionHeight, eyeAzimuth, eyeElevation, eyeMovementGain, screenDistance, screenWidth, screenAspectRatio );

            efferencePriorSigma = 10;
            headingPriorSigma = 10;
            measurementNoiseSigma = 0.0001;

            [params.obsParams, retinalMotionFieldNoiseFree]  = HeadingEstimation.SetupObserver(params.expParams, numberOfRetinalPoints, efferencePriorSigma, headingPriorSigma, measurementNoiseSigma);

            retinalMotionField = retinalMotionFieldNoiseFree;

            [v,w] = HeadingEstimation.FitRegressionModel(retinalMotionField, params);

            table(v, params.expParams.headingVelocity, w, params.expParams.eyeAngularVelocity,  'VariableNames', {'estimatedHeadingVelocity', 'trueHeadingVelocity', 'estimatedEyeAngularVelocity', 'trueEyeAngularVelocity'})
        end

        function RunExampleAnalyticalBiasAndVariance()

        end

        function conditionTable = SwepConditions() 

            headingSpeed    = 10; % m/s
            headingAzimuth  = 00; % deg postiive up

            eyePositionHeight   = 1.5; % m positive up

            numberOfRetinalPoints = 2000;
            screenDistance = 1;
            screenWidth = tand(30)*screenDistance*2;
            screenAspectRatio = 16/9;

            headingPriorSigma = 10;
            measurementNoiseSigma = 1;




            range = 20;
            step = 5;
            az = -range:step:range;
            el = -range/2:step:0;
            gain = [0:1/3:1];
            effCopyNoise = [0.1 10 100 1000];

            % Create grids for all combinations
            [az_grid, el_grid, gain_grid, effCopyNoise_grid] = ndgrid(az, el, gain, effCopyNoise);

            % Combine into a single table
            conditionTable = table(az_grid(:), el_grid(:), gain_grid(:), effCopyNoise_grid(:), ...
                'VariableNames', {'eyeAzimuth', 'eyeElevation', 'eyeMovementGain', 'efferencePriorSigma'});

            
            vx = nan(size(conditionTable.eyeAzimuth));
            vy = nan(size(conditionTable.eyeAzimuth));
            vz = nan(size(conditionTable.eyeAzimuth));
            wx = nan(size(conditionTable.eyeAzimuth));
            wy = nan(size(conditionTable.eyeAzimuth));
            wz = nan(size(conditionTable.eyeAzimuth));

            reps = 100;
            conditionTable = repmat(conditionTable,reps,1);

            height(conditionTable);
            tic
            parfor i=1:height(conditionTable)
                i
                
                params = [];
                eyeAzimuth = conditionTable{i,'eyeAzimuth'};
                eyeElevation = conditionTable{i,'eyeElevation'};
                eyeMovementGain = conditionTable{i,'eyeMovementGain'};
                efferencePriorSigma = conditionTable{i,'efferencePriorSigma'};


                params.expParams = HeadingEstimation.SetUpExperiment(headingSpeed, headingAzimuth, eyePositionHeight, eyeAzimuth, eyeElevation, eyeMovementGain, screenDistance, screenWidth, screenAspectRatio );
                [params.obsParams, retinalMotionFieldNoiseFree]  = HeadingEstimation.SetupObserver(params.expParams, numberOfRetinalPoints, efferencePriorSigma, headingPriorSigma, measurementNoiseSigma);

                retinalMotionField = retinalMotionFieldNoiseFree + randn(size(retinalMotionFieldNoiseFree))*measurementNoiseSigma;

                [v,w] = HeadingEstimation.FitRegressionModel(retinalMotionField, params);

                vx(i) = v(1);
                vy(i) = v(2);
                vz(i) = v(3);
                wx(i) = w(1);
                wy(i) = w(2);
                wz(i) = w(3);

               
            end
toc

conditionTable.vx = vx;
conditionTable.vy = vy;
conditionTable.vz = vz;
conditionTable.wx = wx;
conditionTable.wy = wy;
conditionTable.wz = wz;
conditionTable.headingAzimuth = atan2d(vy,vx);
            a =1
        end


        function [expParams] = SetUpExperiment( headingSpeed, headingAzimuth, eyePositionHeight, eyeAzimuth, eyeElevation, eyeMovementGain, screenDistance, screenWidth, screenAspectRatio )

            expParams.headingSpeed         = headingSpeed;
            expParams.headingAzimuth       = headingAzimuth;
            expParams.headingVelocity      = headingSpeed*[cosd(headingAzimuth) -sind(headingAzimuth) 0]';

            expParams.eyePositionHeight    = eyePositionHeight;
            expParams.eyeAzimuth           = eyeAzimuth;
            expParams.eyeElevation         = eyeElevation;
            expParams.eyeMovementGain      = eyeMovementGain;

            expParams.screenDistance       = screenDistance;
            expParams.screenWidth          = screenWidth;
            expParams.screenAspectRatio    = screenAspectRatio;


            % Calculate eye position and eye velocity

            % get the rotation matrix of the eye (doing a listings rotation so there is
            % no false torsion)
            direction = atan2(eyeElevation, -eyeAzimuth);
            eccentricity = deg2rad(sqrt(eyeElevation.^2 + eyeAzimuth.^2));
            expParams.eyePositionRotMat = Geometry3D.List2Mat([direction eccentricity 0]);
            % get the eye location in the world
            expParams.eyeLocation = [0, 0, eyePositionHeight];

            % calculate the eye linear velocity by multiplying the velocity vector by
            % the rotation matrix of the eye position
            expParams.eyeLinearVelocity = expParams.eyePositionRotMat'*expParams.headingVelocity;

            % calculate the angular eye velocity assuming the eye is pursuing the point
            % in the motion field caused by the heading velocity (times a gain factor)
            foveaDirection = [1 0 0];
            fovealStimAngularVelocity = Geometry3D.CalculateMotionField(foveaDirection, [0 0 0]', expParams.eyeLinearVelocity, expParams.eyeLocation, expParams.eyePositionRotMat);
            [~, JacAngular] = Geometry3D.CalculateMotionJacobianFields(foveaDirection);

            expParams.eyeAngularVelocity = -expParams.eyeMovementGain*JacAngular'*fovealStimAngularVelocity';

            % todo: think about torsion. This formula and the way Jw was calculated
            % implies pursuit following listing's law with axis in the ZY plane
        end

        function [obsParams, retinalMotionField, retinalDepthField] = SetupObserver(expParams, numberOfPoints, efferencePriorSigma, headingPriorSigma, measurementNoiseSigma)

            obsParams.visualDirections = Geometry3D.SampleVisualDirections(numberOfPoints);

            [retinalMotionField, ~, ~, Jv, Jw, retinalDepthField] = Geometry3D.CalculateMotionField(obsParams.visualDirections, expParams.eyeAngularVelocity, expParams.eyeLinearVelocity, expParams.eyeLocation, expParams.eyePositionRotMat);

            obsParams.visualDirectionsJv = Jv;
            obsParams.visualDirectionsJw = Jw;

            obsParams.visualDirectionsDirectionSigma = ones(numberOfPoints,1);% speed dependency?
            obsParams.visualDirectionsSpeedSigma = ones(numberOfPoints,1); % speed dependency?
            obsParams.visualDirectionsVelocityXSigma = measurementNoiseSigma*ones(numberOfPoints,1); % speed dependency?
            obsParams.visualDirectionsVelocityYSigma = measurementNoiseSigma*ones(numberOfPoints,1); % speed dependency?

            obsParams.efferencePriorSigma = efferencePriorSigma; % mean zero
            obsParams.headingPriorSigma = headingPriorSigma; % mean zero
            obsParams.measurementNoiseSigma = measurementNoiseSigma; % mean zero



            % do the stacking of jacobians and multiplying by the depth
            % here to save time

            JacLinear = obsParams.visualDirectionsJv;
            JacAngular = obsParams.visualDirectionsJw;

            % stack the horizontal and vertical measurements and jacobians to do the
            % regression motion = [I D]*[Jang 0; 0; Jline][w v]
            JacLinearStacked = [ squeeze(JacLinear(1,:,:))';squeeze(JacLinear(2,:,:))'];
            JacAngularStacked = [ squeeze(JacAngular(1,:,:))';squeeze(JacAngular(2,:,:))'];
            Jac = [JacAngularStacked zeros(size(JacAngularStacked)); zeros(size(JacAngularStacked)) JacLinearStacked];

            % We need a factor that depends on depth (in diopters) for the linear field
            % but a factor that is just zero or one depending on the presence or
            % absence of stimuli for the angular field. In this case only the floor
            DepthFieldStacked1 = [retinalDepthField zeros(size(retinalDepthField));zeros(size(retinalDepthField)) retinalDepthField];
            DepthFieldStacked = [diag(diag(DepthFieldStacked1)>0) DepthFieldStacked1];

            obsParams.ImageJacobianStacked = DepthFieldStacked*Jac;

        end

        function [estimatedHeadingVelocityHeadRef,estimatedEyeAngularVelocity] = FitRegressionModel(retinalMotionField, params)

            solver = 'Full';

            switch(solver)
                case 'OnlyHeading'

                    % get an estimate of eye velocity (this could be coming from a mix of
                    % efference copy but also motion estimation). For now let's test what
                    % happens changing the gain and adding some noise
                    estimatedEyeVelocityGain = 1;
                    estimatedEyeVelocityNoise = 0;
                    estimatedEyeAngularVelocity = params.expParams.eyeAngularVelocity*estimatedEyeVelocityGain + randn(size(params.expParams.eyeAngularVelocity))*estimatedEyeVelocityNoise;

                    % estimate teh motion field due to that estimated eye velocity and
                    % substract it from the full retinal field to estimate the motion field due
                    % to only linear motion
                    [estimatedAngularMotionField ] = Geometry3D.CalculateMotionField(params.obsParams.visualDirections, estimatedEyeAngularVelocity, [0 0 0]', params.expParams.eyeLocation, params.expParams.eyePositionRotMat);

                    estimatedLinearMotionField = retinalMotionField - estimatedAngularMotionField;
                    estimatedLinearMotionFieldStacked = [estimatedLinearMotionField(:,1);estimatedLinearMotionField(:,2)];

                    % solve via regression
                    estimatedLinearVelocityEyeRef = obsParams.ImageJacobianStacked \ estimatedLinearMotionFieldStacked;

                case 'Full'

                    motionFieldTotalEyeRefefenceStacked = [retinalMotionField(:,1);retinalMotionField(:,2)];

                    efferencePriorSigma = params.obsParams.efferencePriorSigma;
                    headingPriorSigma = params.obsParams.headingPriorSigma;
                    measurementNoiseSigma = params.obsParams.measurementNoiseSigma;


                    % the prior for heading is zero, the prior for angular velocity comes from
                    % the efference copy
                    MuPrior = [params.expParams.eyeAngularVelocity; 0; 0; 0];
                    SigmaPrior = diag([ efferencePriorSigma*ones(1,3) headingPriorSigma*ones(1,3)]);
                    SigmaLikelihood = diag(measurementNoiseSigma*ones(1,height(motionFieldTotalEyeRefefenceStacked)));
                    % Q: should we rotate the prior from eye reference to head reference? Right
                    % now it's in eye reference
                    % Q: how worse would it get if we use a prior for depth instead of the true
                    % depth?

                    % solve via WLS Ridge regression

                    X = params.obsParams.ImageJacobianStacked; % design matrix
                    y = motionFieldTotalEyeRefefenceStacked + randn(size(motionFieldTotalEyeRefefenceStacked))*measurementNoiseSigma; % add noise to measurements

                    estimatedAngularLinearVelocity = (X'/SigmaLikelihood*X + eye(length(MuPrior))/SigmaPrior) \ (X'/SigmaLikelihood*y + SigmaPrior\MuPrior);

                    estimatedLinearVelocityEyeRef = estimatedAngularLinearVelocity(4:6);
                    estimatedEyeAngularVelocity = estimatedAngularLinearVelocity(1:3); % TODO: head or eye reference axis?
            end

            estimatedHeadingVelocityHeadRef = params.expParams.eyePositionRotMat*estimatedLinearVelocityEyeRef;

            
        end



    end
end