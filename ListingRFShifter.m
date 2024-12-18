function shift = ListingRFShifter(eyeAzimuth, eyeElevation, rfAzimuth, rfElevation, listingsPitchDeg, ListingsYawDeg )

    % ListingRFShifter(1) call like this to get a grid plot

    % shift is in meters for a screen at 1 meter


    % We aim to calculate the shift in pixels on a screen between two
    % different scenarios
    % 1 - The vector in screen coordinates between the eye position and the
    % RF position is the same for all eye positions. 
    % 2 - The exact vector in screen coordinates between the eye position
    % and the RF given a rotation of the eye according to Listing's law. 


    if (nargin ~= 1)

        if ( nargin == 0)
            eyeAzimuth = -10;
            eyeElevation = 5;

            rfAzimuth = 20;
            rfElevation = -10;

            listingsPitchDeg = 0;
            ListingsYawDeg = 0;
        end

        % TODO: need to be more specific on what coordinate frame elevation and
        % azimuth are

        % Get the rotation matrix of the gaze direction according to listing's
        % plane

        eyePositionRotMat = Geometry3D.List2Mat([atan2(eyeElevation, -eyeAzimuth) deg2rad(sqrt(eyeElevation.^2 + eyeAzimuth.^2)) 0]);

        % Get the rotation matrix of the RF direction according to listing's
        % plane
        rfPositionRotMat = Geometry3D.List2Mat([atan2(rfElevation, -rfAzimuth) deg2rad(sqrt(rfElevation.^2 + rfAzimuth.^2)) 0]);


        rfVisualDirectionEyeCoords = rfPositionRotMat(:,1);
        rfVisualDirectionHeadCoords = eyePositionRotMat*rfVisualDirectionEyeCoords;

        rfVisualDirectionScreenCoords2 = rfVisualDirectionHeadCoords(2:3)./rfVisualDirectionHeadCoords(1);


        eyeGazeDirection = eyePositionRotMat(:,1);
        eyeGazeDirectionScreenCoords =  eyeGazeDirection(2:3)./eyeGazeDirection(1);
        rfVisualDirectionScreenCoords1 = eyeGazeDirectionScreenCoords + rfVisualDirectionEyeCoords(2:3)./rfVisualDirectionEyeCoords(1);

        shift = rfVisualDirectionScreenCoords1 - rfVisualDirectionScreenCoords2;
    else
        Xeye = [-20:10:20];
        Yeye = [-15:7.5:15];
        Xrf = [-15:3:15];
        Yrf = [-10:2:10];
        [XeyeGrid,YeyeGrid] = meshgrid(Xeye, Yeye);
        [XrfGrid,YrfGrid] = meshgrid(Xrf, Yrf);

        shifts = zeros(length(Xeye), length(Yeye), length(Xrf), length(Yrf), 2);

        for i1=1:length(Xeye)
            for i2=1:length(Yeye)
                for i3=1:length(Xrf)
                    for i4=1:length(Yrf)
                        shifts(i1,i2,i3,i4,:) = ListingRFShifter(Xeye(i1), Yeye(i2), Xrf(i3), Yrf(i4), 0, 0 );
                    end
                end
            end
        end

        shiftsdeg = atand(shifts);

        figure
        tiledlayout(length(Xeye), length(Yeye))
        for i1=1:length(Xeye)
            for i2=1:length(Yeye)
                nexttile
                sx = squeeze(shiftsdeg(i1,i2,:,:,1));
                sy = squeeze(shiftsdeg(i1,i2,:,:,2));
                quiver(XrfGrid(:),YrfGrid(:), sx(:),sy(:), 'off');
                title(sprintf('Eye position (%0.1f,%0.1f)',Xeye(i1), Yeye(i2)))
                set(gca,'xlim',[-20 20],'ylim',[-15 15])
            end
        end
        xlabel('Rf x position (deg)')
        ylabel('Rf y position (deg)')

    end

end