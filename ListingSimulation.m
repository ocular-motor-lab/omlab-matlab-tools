function R = ListingSimulation(eyeAzimuth, eyeElevation)
% This function calculates Rotation matrix that moves the eye from primary
% position to a position where the gaze direction has a particular azimuth
% and elevation while following listing's law.
%
% Azimuth and elevation are defined according to Fick coordinates. That is,
% a gymbal system where horizontal rotation comes first and vertical
% rotation second. 
% IMPORTANTLY we only used them to define the direction of gaze. Not the
% entire 3D position of the eye.
%
%
% We define a right handed coordinate system centered in the eye
%       X is positive forward
%       system)
%       Y is positive left
%       Z is positive up
%
% Then according to the right hand rule rotations are:
%       Horizontal positive towards the left
%       Vertical positive down
%       Torsional positive top to the right
%
% Other coordinate systems may need to be defined unambiguosly too
%
%   Screen coordinate system 
%           Usually 2D cartesian system
%   Camera coordinate system
%           Cartesian coordinate system for rectangular sensor cameras
%           Fick or helmholtz for laser scanning systems
%


% These two rotation matrices are the Fick rotation matrix that produce
% that azimuth and elevation, this rotations will have false torsion. We
% will ignore that. At this point we only care about the direction of gaze.
% That will correspond with the first column of the rotation matrix. 
% 
% for example a rotation of +10 az and +10 elevation
% means to rotate 10 deg to the left and 10 deg down. So the direction of
% goes from [1 0 0] to [0.96 0.17 -0.17]
%
RF = Fick2RotMat(deg2rad([eyeAzimuth, eyeElevation, 0]));

%
% Next we calculate the listing rotation matrix. Here we assume Listing's
% plane is the YZ plane
RL = DirectionOfGazeToListingRotMat(RF(1,1), RF(2,1), RF(3,1));


% Now, this rotation matrices indicate the direcition of gaze (first
% column). The direction of the horizontal meridian of the retina (second
% column). The direction of the vertical meridian of the retina (third
% column). Another way of saying that is that tha plane formed by the
% columns 1 and 2 corresponds with the plane that contains the horizontal
% meridian and the direction of gaze, and columns 1 and 2 define a plane
% that contains the vertical meridian and the direciton of gaze. 

figure
quiver(RL(2,1),(RL(3,1)), RL(2,2),RL(3,2) )
hold
quiver(RL(2,1),(RL(3,1)), RL(2,3),RL(3,3) )
set(gca,'xlim',[-1 1],'ylim',[-1 1])
set(gca, 'XDir', 'reverse');
xlabel('Y axis')
ylabel('Z axis')
end


% Fick to rotation matrix
function M = Fick2RotMat(HVT)

H = HVT(1);
V = HVT(2);
T = HVT(3);

M = RotZ(H)*RotY(V)*RotX(T);
end



% HORIZONTAL ROTATION right handed
function M = RotZ(theta)
M = [   cos(theta)  -sin(theta)     0;
        sin(theta)  cos(theta)      0;
        0           0               1];
end

% VERTICAL ROTATION right handed
function M = RotY(phi)
M = [   cos(phi)  0               sin(phi);
        0           1               0;
        -sin(phi) 0               cos(phi)];
end

% TORSIONAL ROTATION right handed
function M = RotX(psi)
M = [   1           0               0;
        0           cos(psi)      -sin(psi);
        0           sin(psi)      cos(psi)];
end

function M = DirectionOfGazeToListingRotMat(x,y,z)

%
% To do this we first calculate a sort of  polar coordinates of the
% direction of gaze. What is the eccentricity and what is the direction. 
%
% In the example (+10, +10) this corresponds to an eccentricity of 14 deg
% and an angle of 45. This angle will be zero for horizontal movement to
% the left, and will be 90 deg for a vertical movement up. 

ecc = acos(x);
angle = atan2(-z,y);
M = List2Mat([angle ecc 0]);
end


function M = List2Mat(ADT)
% the input is angle, eccentricity and torsion
% sequency of rotations is as follows
% 1- rotate the coordinate system according the angle (
% then rotate an eccentricity ammount arond the rotated Z axis
% finally undo the intial torsion and add the actual torsion
A = ADT(1);
D = ADT(2);
T = ADT(3);
M = RotX(T-A)*RotZ(D)*RotX(A);
end