%% SO3 attractor simulation
%% Jorge Otero-Millan 3/10/2024
close all
N = 4; % 2 ring, 3 sphere, 4 quaternions (SO3);

% time parameters of the simulation
dt = 0.001;
t = (0:dt:20)';

% initial conditions of the attractor
x0 = 2*[1 0 0 0 1 0 0 0]'; % point outside the attractor to show initial drift towards it

% velocity input
v = zeros(length(t),3); 
v(t > 1 & t <2, 1) = 1; % angular velocity around x
v(t > 3 & t <4, 2) = 5; % angular velocity around y
v(t > 5 & t <9, 3) = -1; % angular velocity around z

% fixed point the system drift towards in the absence of input. 
xf = 0*[1 0 0 0]'; % multiply by  gain of the drift towards the fixed point.

% position inputs (we have two to test bayesian integrator, together with
% the fixed point as a prior)
p = zeros(length(t),4); 
Gp1 = 0;
Gp2 = 0;
p(t > 16 & t <20, :) = ... 
    Gp1*repmat(eul2quat(deg2rad([30,0,0])), sum(t > 16 & t <20),1)+ ...
    Gp2*repmat(eul2quat(deg2rad([20,0,0])), sum(t > 16 & t <20),1);

% so3 attractor matrix definition, conic section extended to be the
% hypersphere containing SO3. x'Ax =0 if x1^2 + x2^2 + x3^2 + x4^2 - 1 = 0
A = eye(5);
A(end,end)=-1;

% base of the subspace that contains the attractor
S = [eye(4);eye(4)];


% A tensor represents the group action inducing moving towards the
% attractor (first row) and along the atractor (next rows).
T = cat(3, ...
    [1 0 0 0 0;    0 1 0 0 0;  0 0 1 0 0;  0 0 0 1 0;  0 0 0 0 0], ...
    [0 -1 0 0 0;   1 0 0 0 0;  0 0 0 -1 0; 0 0 1 0 0;  0 0 0 0 0], ...
    [0 0 -1 0 0;   0 0 0 1 0;  1 0 0 0 0;  0 -1 0 0 0; 0 0 0 0 0], ...
    [0 0 0 -1 0;   0 0 -1 0 0; 0 1 0 0 0;  1 0 0 0 0;  0 0 0 0 0]);

[t, xout] = ode45(@(ti,xi)odeAttractor4D(...
    ti, xi, ...
    interp1(t,v,ti)', interp1(t,p,ti)', ...
    A,T, xf), ...
    t, x0);

PlotRun(t,v,xout);

function xd = odeAttractor4D( ti, xi, vi, pi, A, T, S, xf)

    n = length(xi); % dimensions space
    np = size(S,2); % dimensions of the subspace containing the manifold
    P = inv(S'*S)*S'; % projection matrix to the plane
    
    P = [P zeros(np,1); zeros(1,n) 1]; % homogenous projection matrix to allow translation of the plane away from zero
    S = [S zeros(n,1); zeros(1,np) 1];


    x   = [xi;1]; % make it homogeneous
    xf  = [xf;1]; 
    p   = [pi;1];
    v   = vi;

    xx = x;
    x = P*x;
    

    % Tensor product to rotate an orthogonal frame and align it with the
    % surface of the attractor.
    % Dim1 is the direction of movement towards the attractor
    % Dims 1..n are the directions of movement along the attractor dim1
    % driven by the input 
    M = zeros(length(x),size(T,3));
    for i=1:size(T,3)
        M(:,i) = T(:,:,i) * A*x;
    end

    % Scale the directions of motion by the velocity
    xd = M * [ - (4*x'*A*x);    ...    % velocity towards the attractor
              v ]               ...    % velocity along the attractor scaled by velocity input
         + (xf - x)             ...    % velocity towards the fixed point (prior)
         + (p - x);             ...    % velocity towards the position input (likelihood)

    % remove the homogenous component
    xd = S*xd + 10*(S*x-xx); % added velocity towards the manifold subspace
    xd = xd(1:end-1);
end

function PlotRun(t,x,xout)
    figure
    subplot(3,3,1);
    c = linspace(1,255,length(xout));
    scatter(xout(:,1),xout(:,2),[],c); set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
    xlabel( 'unit 1'), ylabel( 'unit 2')
    title('SO3 attractor')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(3,3,4);
    scatter(xout(:,2),xout(:,3),[],c);  set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
    xlabel( 'unit 2'), ylabel( 'unit 3')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(3,3,7);
    scatter(xout(:,3),xout(:,4),[],c);  set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
    xlabel( 'unit 3'), ylabel( 'unit 4')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    

    subplot(3,3,[2 3]);
    plot(t,x,'linewidth',2);hold
    plot(t,xout);
    legend({'Input velocity x', 'Input velocity y', 'Input velocity z', 'unit 1', 'unit 2', 'unit 3', 'unit 4'})
    xlabel('Time')

    subplot(3,3,[5 6]);
    angles = rad2deg(quat2eul(xout(:,1:4),'XYZ'));
    plot(t,angles);
    set(gca,'clim', [min(min(xout(t>1,:))), max(max(xout(t>1,:)))]) % make the clim ignore the first 100 timepoints
    xlabel('Time')
    ylabel('angle(deg)')
    % set(gca,'xticklabel',[],'yticklabel',[])


    subplot(3,3,[8 9]);
    imagesc(xout')
    set(gca,'clim', [min(min(xout(t>1,:))), max(max(xout(t>1,:)))]) % make the clim ignore the first 100 timepoints
    xlabel('Time')
    ylabel('units')
    set(gca,'xticklabel',[],'yticklabel',[])
end
