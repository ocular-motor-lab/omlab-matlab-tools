%% attractor simulation
clear params;
close all

% time parameters of the simulation
tstart = 0;
tend = 20;
dt = 0.001;
trange = (tstart:dt:tend)';

% initial conditions of the attractor
y0 = [0 0]';

% input velocity
x = zeros(size(trange));
x(trange > 1 & trange <2) = 10;
x(trange > 3 & trange <4) = 20;
x(trange > 5 & trange <9) = -10;
x = x+ randn(size(x))*.1;

% ring attractor

A = 1;
B = 0;
C = 1;
D = -2;
E = -2;
F = 1;
params.Aq  = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
params.FP = [1 1] - [cosd(45) sind(45)];


[t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x,params), trange, y0);

figure
subplot(3,3,1);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
xlabel( 'Internal unit 1'), ylabel( 'Internal unit 2')
title('Ring attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[2 3]);
plot(t,[x e]);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')


% parabollic attractor

A = 0;
B = -1;
C = 1;
D = -1;
E = -1;
F = 1;
params.Aq = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
params.FP = -[0.5 0.5];

[t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x/10,params), trange, y0); % velocity made slower to avoid numerical problems

subplot(3,3,4);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Parabollic attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[5 6]);
plot(t,[x e]);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')

% % line attractor
A = 0;
B = 0;
C = 0;
D = 1;
E = 1;
F = -(D+E);
params.Aq  =  [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
params.FP = [1 1];

[t, e] = odeSimple(@(t,y)ode2DAttractor(t,y,trange,x,params), trange, y0);

subplot(3,3,7);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Line attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[8 9]);
plot(t,[x e]);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')


function yd = ode2DAttractor(t,y,xt,x,params)
    x = interp1(xt,x,t); % Interpolate the data set (gt,g) at time t
    
    y = [y;1]; % make it homogeneous

    % x is the input
    % y is the state/output of the atractor

    W = params.Aq;
    Ortho1 = [0 -1 0;1 0 0; 0 0 0]; %rotates Wy to get the orthogonal direction along the attractor
    Gleak = 0.2; % leak speed

    fp = params.FP;
    leak = y(2)*fp(1)-y(1)*fp(2);
    
    yd = (x-Gleak*leak)*Ortho1*W*y   -   (4*y'*W*y)*W*y;

    yd = yd(1:2);
    
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


function [t, yout] = odeSimple(fun, trange, y0)
if ( length(trange) == 2 )
    trange = linspace(trange(1),trange(2),1000);
end
yout= zeros(length(y0),length(trange));
t = trange;

yout(:,1) = y0;
for i=2:length(t)
    yout(:,i) = yout(:,i-1) + fun(t(i), yout(:,i-1))*(t(i)-t(i-1)) ;
end

yout = yout';

end

