%% SO3 attractor simulation
%% Jorge Otero-Millan 3/10/2024
close all

% time parameters of the simulation
dt = 0.001;
t = (0:dt:20)';

% initial conditions of the attractor
y0 = [1;1;0;0]; % point outside the attractor to show initial drift towards it

% input velocity
x = zeros(length(t),3);
x(t > 1 & t <2, 1) = 1; % angular velocity around x
x(t > 3 & t <4, 2) = 2; % angular velocity around y
x(t > 5 & t <9, 3) = -1; % angular velocity around z

% so3 attractor matrix definition, conic section extended to be the
% hypersphere containing SO3. x'Ax =0 if x1^2 + x2^2 + x3^2 + x4^2 - 1 = 0
A = eye(5);
A(end,end)=-1;

% A tensor represents the group action inducing moving towards the
% attractor (first row) and along the atractor (next rows).
T = cat(3, ...
    [1 0 0 0 0;    0 1 0 0 0;  0 0 1 0 0;  0 0 0 1 0;  0 0 0 0 0], ...
    [0 -1 0 0 0;   1 0 0 0 0;  0 0 0 -1 0; 0 0 1 0 0;  0 0 0 0 0], ...
    [0 0 -1 0 0;   0 0 0 1 0;  1 0 0 0 0;  0 -1 0 0 0; 0 0 0 0 0], ...
    [0 0 0 -1 0;   0 0 -1 0 0; 0 1 0 0 0;  1 0 0 0 0;  0 0 0 0 0]);

[t, yout] = ode45(@(ti,yi)odeAttractor4D(yi,interp1(t,x,ti),A,T), t, y0);

PlotRun(t,x,yout);

function yd = odeAttractor4D( y, x, A, T)
    
    y = [y;1]; % make it homogeneous

    yd = - (4*y'*A*y) * T(:,:,1) * A*y ... % movement towards the attractor
               + x(1) * T(:,:,2) * A*y ... % movement along the attractor dim1 driven by the input
               + x(2) * T(:,:,3) * A*y ... % movement along the attractor dim2 driven by the input
               + x(3) * T(:,:,4) * A*y;    % movement along the attractor dim3 driven by the input

    yd = yd(1:end-1); % remove the homogenous component
end

function PlotRun(t,x,yout)
    figure
    subplot(3,3,1);
    c = linspace(1,255,length(yout));
    scatter(yout(:,1),yout(:,2),[],c); set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
    xlabel( 'unit 1'), ylabel( 'unit 2')
    title('SO3 attractor')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(3,3,4);
    scatter(yout(:,2),yout(:,3),[],c);  set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
    xlabel( 'unit 2'), ylabel( 'unit 3')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(3,3,7);
    scatter(yout(:,3),yout(:,4),[],c);  set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
    xlabel( 'unit 3'), ylabel( 'unit 4')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(3,3,[2 3]);
    plot(t,x,'linewidth',2);hold
    plot(t,yout);
    legend({'Input velocity x', 'Input velocity y', 'Input velocity z', 'unit 1', 'unit 2', 'unit 3', 'unit 4'})
    xlabel('Time')
    subplot(3,3,[5 6 8 9]);
    imagesc(yout')
    set(gca,'clim', [min(min(yout(t>1,:))), max(max(yout(t>1,:)))]) % make the clim ignore the first 100 timepoints
    xlabel('Time')
    ylabel('units')
    set(gca,'xticklabel',[],'yticklabel',[])
end
