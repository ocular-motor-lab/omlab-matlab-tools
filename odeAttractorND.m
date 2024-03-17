function dx = odeAttractorND( ti, xi, wi, A, T, S)
% ODEATTRACTORND  differential equation for an attractor network of m
% dimensions embedded in an n dimensional space with p inputs. 
%
%   dx = odeAttractorND( ti, xi, wi, A, T, S) 
%
%   outputs:
%       dx: derivative of the state x at time ti
%
%   inputs: 
%       ti: (1 x 1) current time ti (very useful for debugging)
%       xi: (n x 1) state at time ti 
%       wi: (q x 1) velocity input at time ti 
%       A:  (m+1 x m+1) matrix defining the quadratic ecuation that specifis
%           the attractor shape. The attractor is the set of points x that
%           meet the condition x'*A*x=0
%
%           For example:
%           for 1D rotations:
%               A = [1 0 0; ...
%                    0 1 0; ...
%                    0 0 -1];
%           for 2D rotations:
%               A = [1 0 0 0; ...
%                    0 1 0 0; ...
%                    0 0 1 0; ...
%                    0 0 0 -1];
%           for 3D rotations:
%               A = [1 0 0 0 0; ...
%                    0 1 0 0 0; ...
%                    0 0 1 0 0; ...
%                    0 0 0 1 0; ...
%                    0 0 0 0 -1];
%
%       T: (m+1 x m+1 x q+1) Tensor encoding a base of the skew symmetric
%           matrices of size m. Each of the q matrices roates the direction
%           towards the attractor so we get one direction q directions
%           along the tangent of the attractor. 
%
%           For example:
%           for 1D rotations:
%               T = cat(3, ...
%                   [0 -1 0 ; 1 0 0 ; 0 0 0;]); % 3x3x2
%
%           for 2D rotations:
%               T = cat(3, ...
%                   [0 -1 0 0;  1 0 0 0;  0 0 0 0;  0 0 0 0], ...
%                   [0 0 1 0;   0 0 0 0;  -1 0 0 0; 0 0 0 0]); % 4x4x3
%
%           for 3D rotations:
%               T = cat(3, ...
%                   [0 -1 0 0 0;   1 0 0 0 0;  0 0 0 -1 0; 0 0 1 0 0;  0 0 0 0 0], ...
%                   [0 0 -1 0 0;   0 0 0 1 0;  1 0 0 0 0;  0 -1 0 0 0; 0 0 0 0 0], ...
%                   [0 0 0 -1 0;   0 0 -1 0 0; 0 1 0 0 0;  1 0 0 0 0;  0 0 0 0 0]); % 5x5x3
%
%       S: (n x m) Matrix to project x into a subspace that contains the
%           attractor. Each column is a base vector defining the subspace.
%
%
%   Jorge Otero-Millan
%   Ocular-Motor lab, UC Berkeley
%   3/16/2024
%
%
    x = [xi;1]; % make the state homogeneous
    w = wi;     % just change the name for consistency
    [n, m] = size(S); % dimensions space and subspace
    
    % Calculate the projection matrix subspace S that contains the
    % attractor. Added column and row for the homogenous component 
    P = [(S'*S)\S' zeros(m,1); zeros(1,n) 1];
    S = [S zeros(n,1); zeros(1,m) 1];

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

    dx = dx(1:end-1); % remove the homogenous component
end
