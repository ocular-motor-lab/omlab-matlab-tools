function xd = odeAttractorND( ti, xi, wi, pi, A, T, S, xf)
% ODEATTRACTORND  differential equation for an attractor network of m
% dimensions embedded in an n dimensional space with p inputs. 
%   xd = odeAttractorND( ti, xi, vi, pi, A, T, S, xf) 
%
%   outputs:
%       xd: derivative of the state x at time ti
%
%   inputs: 
%       ti: (1 x 1) current time ti
%       xi: (n x 1) state at time ti 
%       wi: (q x 1) velocity input at time ti 
%       pi: (n x 1) position input at time ti 
%       A:  (m+1 x m+1) matrix defining the quadratic ecuation that specifis
%           the attractor shape. The attractor is the set of points x that
%           meet the condition x'*A*x=0
%       T: (m+1 x m+1 x q+1) Tensor encoding the base of the skew sydirections of
%           motion from x relative to the direction orthogonal to the
%           attractor. Each of the p+1 matrices roates the direction
%           towards the attractor so we get one direction typically towards
%           the attractor and p directions along the attractor.
%       S: (n x m) Matrix to project x into a subspace that contains the
%           attractor. Each column is a base vector defining the subspace.
%       xf: (n x 1) Fixed point (set point or null point) on the space. The
%           state will drift towards it in the absence of inputs.
%

    [n, m] = size(S); % dimensions space and subspace
    q = size(T,3)-1; % number of input dimensions
    
    % P = inv(S'*S)*S'; % projection matrix to the plane if S not
    % orthonormal
    P = [S' zeros(m,1); zeros(1,n) 1]; % homogenous projection matrix to allow translation of the plane away from zero
    S = [S zeros(n,1); zeros(1,m) 1];

    x  = [xi;1]; % make it homogeneous
    xf = [xf;1]; 
    p  = [pi;1];
    w  = wi;
    
    % vector pointing towards the attractor
    v = A*P*x; 

    % Tensor product to rotate an orthogonal frame and align it with the
    % surface of the attractor.
    % Dim1 is the direction of movement towards the attractor
    % Dims 1..n are the directions of movement along the attractor dim1
    % driven by the input 
    C = squeeze(pagemtimes(T, v));

    % Scale the directions of motion by the velocity
    xd = C * w;                      % velocity along the attractor scaled by velocity input
    xd = xd - 10*(4*x'*P'*A*P*x*v); ... % velocity towards the attractor
    
%     xd = xd + norm(xf)*(P*xf - P*x)  ...    % velocity towards the fixed point (prior)
%             + norm(p)*(P*p - P*x);   ...    % velocity towards the position input (likelihood)
%             
    xd = xd + norm(xf)*C*diag([0 ones(1,q)])*C'*(P*xf - P*x); ... % velocity towards the fixed point (prior)
    xd = xd + norm(p)*C*diag([0 ones(1,q)])*C'*(P*p  - P*x); ...  % velocity towards the position input (likelihood)

    xd = S*xd + 10*(S*P*x-x); % added velocity towards the manifold subspace
    xd = xd(1:end-1); % remove the homogenous component
end

% IDEA for the drifts

% if  Ax and (xf-x) are parallel it means that we don't need to move