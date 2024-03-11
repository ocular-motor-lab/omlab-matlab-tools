%% SO3 attractor simulation
%% Jorge Otero-Millan 3/10/2024


clear params;
close all

% time parameters of the simulation
tstart = 0;
tend = 20;
dt = 0.001;
trange = (tstart:dt:tend)';

% initial conditions of the attractor
y0 = [1 0 0 0]';

% input velocity
x = zeros(size(trange));
x(trange > 1 & trange <2) = 1; % angular velocity around x
y = zeros(size(trange));
y(trange > 3 & trange <4) = 2; % angular velocity around y
z = zeros(size(trange));
z(trange > 5 & trange <9) = -1; % angular velocity around y
x =[x y z];

% so3 attractor
params.Aq = eye(5);
params.Aq(end,end)=-1; % conic section extended to be the hypersphere containing SO3


[t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x, params), trange, y0);

figure
subplot(3,3,1);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2)
xlabel( 'unit 1'), ylabel( 'unit 2')
title('SO3 attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])

subplot(3,3,4);
plot(e(:,2),e(:,3)); set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2)
xlabel( 'unit 2'), ylabel( 'unit 3')
set(gca,'PlotBoxAspectRatio',[1 1 1])

subplot(3,3,7);
plot(e(:,3),e(:,4)); set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2)
xlabel( 'unit 3'), ylabel( 'unit 4')
set(gca,'PlotBoxAspectRatio',[1 1 1])


subplot(3,3,[2 3]);
plot(t,x,'linewidth',2);hold
plot(t,e);
legend({'Input velocity x', 'Input velocity y', 'Input velocity z', 'unit 1', 'unit 2', 'unit 3', 'unit 4'})
xlabel('Time')
subplot(3,3,[5 6 8 9]);
imagesc(e')
set(gca,'clim', [min(min(e(trange>1,:))), max(max(e(trange>1,:)))]) % make the clim ignore the first 100 timepoints
xlabel('Time')
ylabel('units')
set(gca,'xticklabel',[],'yticklabel',[])


function yd = ode2DAttractor(t,y,xt,x,params)
    x = interp1(xt,x,t); % Interpolate the data set (gt,g) at time t
    
    y = [y;1]; % make it homogeneous

    % x is the input
    % y is the state/output of the atractor

    W = params.Aq;
    Ortho1 = [0 -1 0 0 0;1 0 0 0 0; 0 0 0 -1 0; 0 0 1 0 0; 0 0 0 0 0]; %rotates Wy to get the orthogonal directions along the attractor
    Ortho2 = [0 0 -1 0 0;0 0 0 1 0;1 0 0 0 0; 0 -1 0 0 0; 0 0 0 0 0]; 
    Ortho3 = [0 0 0 -1 0;0 0 -1 0 0;0 1 0 0 0; 1 0 0 0 0; 0 0 0 0 0]; 

    yd = x(:,1)*Ortho1*W*y +x(:,2)*Ortho2*W*y +x(:,3)*Ortho3*W*y   -   (4*y'*W*y)*W*y;

    yd = yd(1:end-1);
    
    % the first term is the movement along the attractor. It can be either
    % driven by the input or by a leak towards a point in the attractor.

    % right now the leak component is not ideal and does not extend to 3D
    % it is sort of a cross product between the current position along the
    % attractor and the set point. So you get the velocity that you need to
    % drift by depending on how far you are. In this case I am assuming
    % the set point is along the diagonal. But that is something that
    % should be done via weights that need to learned too. 

    % the second term is the one that attracts to the conic curve
    % which can be a line, a circle or whatever. 

    % the first term needs to learn the cross product or quaternion product
    % (when we move to more dimensions)
    % for more dimensions we need to calculate the cross product (somehow
    % corrected by manifold shape) between the fixeed point and the current
    % point. That will give the leak velocity which can be added in all
    % components to x. Then x needs to be a applied along the orthogonal
    % dimensions to the Aq*yb which is a vector towards the manifold. So
    % those orthogonal components are the components along the manifold. 

    % I stil think there may be a way to combine the two components with a
    % single tensor. The tensor will define how much to movement a long
    % each of the components. 1 towards the attractor and n-1 along the
    % attractor + 1 along the input regarless of attractor

end

