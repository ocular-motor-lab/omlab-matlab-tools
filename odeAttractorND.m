function xd = odeAttractorND( ti, xi, wi, A, T, S)
% ODEATTRACTORND  differential equation for an attractor network of m
% dimensions embedded in an n dimensional space with p inputs. 
%   xd = odeAttractorND( ti, xi, wi, A, T, S) 
%
%   outputs:
%       xd: derivative of the state x at time ti
%
%   inputs: 
%       ti: (1 x 1) current time ti
%       xi: (n x 1) state at time ti 
%       wi: (q x 1) velocity input at time ti 
%       A:  (m+1 x m+1) matrix defining the quadratic ecuation that specifis
%           the attractor shape. The attractor is the set of points x that
%           meet the condition x'*A*x=0
%
%           For example,
%           for 1D rotations:
%               A = [1 0 0; ...
%                    0 1 0; ...
%                    0 0 -1];
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
%           For example,
%           for 1D rotations:
%               T = cat(3, ...
%                   [0 -1 0 ;   1 0 0 ;  0 0 0;]);
%
%           for 3D rotations:
%               T = cat(3, ...
%                   [0 -1 0 0 0;   1 0 0 0 0;  0 0 0 -1 0; 0 0 1 0 0;  0 0 0 0 0], ...
%                   [0 0 -1 0 0;   0 0 0 1 0;  1 0 0 0 0;  0 -1 0 0 0; 0 0 0 0 0], ...
%                   [0 0 0 -1 0;   0 0 -1 0 0; 0 1 0 0 0;  1 0 0 0 0;  0 0 0 0 0]);
%
%       S: (n x m) Matrix to project x into a subspace that contains the
%           attractor. Each column is a base vector defining the subspace.
%

    [n, m] = size(S); % dimensions space and subspace
    
    % Calculate the projection matrix to the plane S. Make it homogenous
    % projection matrix to allow translation of the plane away from zero 
    P = [(S'*S)\S' zeros(m,1); zeros(1,n) 1];
    S = [S zeros(n,1); zeros(1,m) 1];

    x  = [xi;1]; % make the state homogeneous
    w  = wi;

    % Tensor product to scale the base of the tangent space. Then multiply by That is go
    % from the vector APx that points towards the attractor to a 
    % an orthogonal bases that is tangent to the attractor so it allows for
    % movement along it's surface 
if  (ti>8)
a=1;
end


%     tensorprod(T,A*P*x,2,1)*w;
%     tensorprod(T,w,3,1)* A*P*x ;
    xd = tensorprod(T,w,3,1) * A*P*x ... % velocity along the attractor scaled by velocity input
            - (4*x'*P'*A*P*x*A*P*x); ... % velocity towards the attractor

    xd = S*xd + (S*P*x-x); % added velocity towards the manifold subspace 
    xd = xd(1:end-1); % remove the homogenous component
end
