%% SO3 attractor simulation
%% Jorge Otero-Millan 3/10/2024
close all

% time parameters of the simulation
dt = 0.001;
t = (0:dt:20)';

% initial conditions of the attractor
x0 = [1 1 0 0]'; % point outside the attractor to show initial drift towards it

% input velocity
v = zeros(length(t),3); 
v(t > 1 & t <2, 1) = 1; % angular velocity around x
v(t > 3 & t <4, 2) = 5; % angular velocity around y
v(t > 5 & t <9, 3) = -1; % angular velocity around z

% position inputs
p1 = zeros(length(t),4); 
p1(t > 16 & t <20, :) = repmat(eul2quat(deg2rad([30,0,0])), sum(t > 16 & t <20),1);
p2 = zeros(length(t),4); 
p2(t > 16 & t <20, :) = 10*repmat(eul2quat(deg2rad([20,0,0])), sum(t > 16 & t <20),1);

% so3 attractor matrix definition, conic section extended to be the
% hypersphere containing SO3. x'Ax =0 if x1^2 + x2^2 + x3^2 + x4^2 - 1 = 0
A = eye(5);
A(end,end)=-1;

% fixed point the system drift towards in the absence of input. 
xf = [1 0 0 0]';
Gf = 1; % gain of the drift towards the fixed point.

% A tensor represents the group action inducing moving towards the
% attractor (first row) and along the atractor (next rows).
T = cat(3, ...
    [1 0 0 0 0;    0 1 0 0 0;  0 0 1 0 0;  0 0 0 1 0;  0 0 0 0 0], ...
    [0 -1 0 0 0;   1 0 0 0 0;  0 0 0 -1 0; 0 0 1 0 0;  0 0 0 0 0], ...
    [0 0 -1 0 0;   0 0 0 1 0;  1 0 0 0 0;  0 -1 0 0 0; 0 0 0 0 0], ...
    [0 0 0 -1 0;   0 0 -1 0 0; 0 1 0 0 0;  1 0 0 0 0;  0 0 0 0 0]);

[t, xout] = ode45(@(ti,xi)odeAttractor4D(...
    ti, xi, ...
    interp1(t,v,ti)', interp1(t,p1,ti)', interp1(t,p2,ti)', ...
    A,T, xf, Gf), ...
    t, x0);

PlotRun(t,v,xout);

function xd = odeAttractor4D( t, x, v, p1, p2, A, T, xf, Gf)
    if ( t > 16)
        a=1;
    end
    
    x = [x;1]; % make it homogeneous
    
    M = [T(:,:,1) * A*x, ... % direction of movement towards the attractor
         T(:,:,2) * A*x, ... % direction of movement along the attractor dim1 driven by the input
         T(:,:,3) * A*x, ... % direction of movement along the attractor dim2 driven by the input
         T(:,:,4) * A*x];    % direction of movement along the attractor dim3 driven by the input
    
    % distance between the current point and the fixed point. Ideally this
    % should be the difference in integral along the manifold. But I don't
    % know how to do that. For now, as long as the manifold is a conic this
    % should aproximately work.
    % Then project that distance into M directions to see how much we need
    % to move along each of them. Effectively the dot product
    driftv = Gf*M' * ([xf;1] - x);

    % same thing for each of the position inputs
    driftv = driftv + M' * ([p1;1] - x) + M' * ([p2;1] - x);

    % The total velocity now is the velocity input and the velocity towards
    % the fixed point and the velocity towards 
    totalv = [ - (4*x'*A*x); ... % velocity towards the attractor
              v + driftv(2:end)]; % velocity along the attractor

    % Scale the directions of motion by the velocity
    xd = M * totalv;
    
    % remove the homogenous component
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
    angles = rad2deg(quat2eul(xout,'XYZ'));
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
