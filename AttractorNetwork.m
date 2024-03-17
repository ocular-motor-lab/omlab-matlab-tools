function dx = AttractorNetwork( t, x, w, A, T, S)
% ODEATTRACTORND  differential equation for an attractor network of m
% dimensions embedded in an n dimensional space with m-1 inputs. 
%
%   dx = odeAttractorND( t, x, w, A, T, S) 
%
%   outputs:
%
%       dx: derivative of the state x at time ti
%
%   inputs: 
%
%       t: (1 x 1) current time (very useful for debugging)
%
%       x: (n x 1) state at time t 
%
%       w: (m-1 x 1) velocity input at time t 
%
%       A:  (m+1 x m+1) matrix defining the quadratic ecuation that specifis
%           the attractor shape. The attractor is the set of points x that
%           meet the condition x'*A*x=0
%
%           For example:
%
%           for 1D rotations:
%               A =   1   0   0
%                     0   1   0
%                     0   0  -1
%
%           for 2D rotations:
%               A =   1   0   0   0 
%                     0   1   0   0 
%                     0   0   1   0 
%                     0   0   0  -1  
%
%           for 3D rotations:
%               A =   1   0   0   0   0 
%                     0   1   0   0   0 
%                     0   0   1   0   0 
%                     0   0   0   1   0 
%                     0   0   0   0  -1 
%
%       T: (m x m x m-1) Tensor containing a basis for the m+1 m+1 skew
%           symmetric matrices of size m. Each of the q matrices roates the direction
%           towards the attractor so we get one direction q directions
%           along the tangent of the attractor. 
%
%           For example:
%
%           for 1D rotations:
%                T(:,:,1) =     0    -1
%                               1     0
% 
%           for 2D rotations:
%
%                T(:,:,1) =     0    -1     0
%                               1     0     0
%                               0     0     0
% 
%                T(:,:,2) =     0     0     1
%                               0     0     0
%                              -1     0     0
%
%           for 3D rotations:
%
%                T(:,:,1) =     0    -1     0     0
%                               1     0     0     0
%                               0     0     0    -1
%                               0     0     1     0
% 
%                T(:,:,2) =     0     0    -1     0
%                               0     0     0     1
%                               1     0     0     0
%                               0    -1     0     0
% 
%                T(:,:,3) =     0     0     0    -1
%                               0     0    -1     0
%                               0     1     0     0
%                               1     0     0     0
% 
%       S: (n x m) Matrix to project x into a subspace that contains the
%           attractor. Each column is a base vector defining the subspace.
%           If using n=m just make it the identity.
%
%       Example running of ring attractor:
%           
%            t = (0:0.001:5)';
%            A  = [1 0 0; 0 1 0; 0 0 -1 ];
%            T = [0 -1;   1 0 ; ];
%            S = [1 0; 0 1];
%            w = zeros(size(t));
%            w(t >= 2 & t <4) = deg2rad(45); % 45 deg/s for 2 seconds
%            x0 = [1 0]';
%            [t, xout] = ode45(@(ti,xi) ...
%                    odeAttractorND( ti, xi, interp1(t,w,ti)', A, T, S), ...
%                    t, x0);
%            figure, plot(t,xout);
% 
%
%   Jorge Otero-Millan
%   Ocular-Motor lab, UC Berkeley
%   3/16/2024
%
%
    x = [x;1]; % make the state homogeneous

    [n, m] = size(S); % dimensions space and subspace
    
    % Calculate the projection matrix for subspace S that contains the
    % attractor
    P = (S'*S)\S';
    
    % Pad the matrices to deal with the homogenous component
    P = [ P   zeros(m,1)     ;  zeros(1,n)    1  ];
    S = [ S   zeros(n,1)     ;  zeros(1,m)    1  ];
    T = [ T   zeros(m,1,m-1) ;  zeros(1,m+1,m-1) ];

    %
    % Differential equation dx/dt = f(x, w | A, T, S)
    %
    % First term: (S*tensorprod(T,w,3,1) * A*P*x )
    % Tensor product to scale the base of the tangent space. Then multiply
    % to go from the vector APx that points towards the attractor to a 
    % an orthogonal basis that is tangent to the attractor allowing for
    % movement along its surface depending on the velocity
    %
    % Second term: S*(4*x'*P'*A*P*x*A*P*x) 
    % Velocity towards the attractor following the gradient of the squared
    % quadratic equation.
    %
    % Third term: (S*P*x-x)
    % Added velocity towards the manifold subspace 

    dx = S*tensorprod(T,w,3,1) * A*P*x - S*(4*x'*P'*A*P*x*A*P*x) + (S*P*x-x);

    dx = dx(1:end-1); % remove the homogeneous component
end
