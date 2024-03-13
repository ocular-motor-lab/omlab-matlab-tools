%% attractor simulation
clear params;
% close all

% time parameters of the simulation
dt = 0.001;
t = (0:dt:20)';

% initial conditions of the attractor
x0 = [0 0]';

% input velocity
v = zeros(size(t));
v(t > 1 & t <2) = 10;
v(t > 3 & t <4) = 20;
v(t > 5 & t <9) = -10;
v = v+ randn(size(v))*.1;

% ring attractor
A = 1;
B = 0;
C = 1;
D = -2;
E = -2;
F = 1;
Aq  = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
xf = 0.1*([1 1]' - [cosd(45) sind(45)]');
n=2;
S = eye(n);
p = zeros(length(t),n);
T = cat(3, ...
    [1 0 0 ;    0 1 0 ;  0 0 0;], ...
    [0 -1 0 ;   1 0 0 ;  0 0 0;]);
Gleak = 0.2;

[t, xout] = ode45(@(ti,xi)odeAttractorND(...
    ti, xi, ...
    interp1(t,v,ti)', interp1(t,p,ti)', ...
    Aq,T,S, xf), ...
    t, x0); 

figure
    c = linspace(1,255,length(xout));
subplot(3,3,1);
scatter(xout(:,1),xout(:,2),[],c); set(gca,'xlim',[-1 3]*1.2,'ylim',[-1 3]*1.2), colormap(gca,'jet')
xlabel( 'Internal unit 1'), ylabel( 'Internal unit 2')
title('Ring attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[2 3]);
plot(t,[v xout]);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')


% parabollic attractor

A = 0;
B = -1;
C = 1;
D = -1;
E = -1;
F = 1;
Aq = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
xf = -0.1*[0.5 0.5]';

[t, xout] = ode45(@(ti,xi)odeAttractorND(...
    ti, xi, ...
    interp1(t,v/5,ti)', interp1(t,p,ti)', ...
    Aq,T,S, xf), ...
    t, x0); % velocity made slower to avoid numerical problems

subplot(3,3,4);
scatter(xout(:,1),xout(:,2),[],c); set(gca,'xlim',[-1 3]*1.2,'ylim',[-1 3]*1.2), colormap(gca,'jet')
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Parabollic attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[5 6]);
plot(t,[v xout]);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')

% % line attractor
A = 0;
B = 0;
C = 0;
D = 1;
E = 1;
F = -(D+E);
Aq  =  [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
xf = 0.1*[1 1]';

S = eye(2);

[t, xout] = ode45(@(ti,xi)odeAttractorND(...
    ti, xi, ...
    interp1(t,v,ti)', interp1(t,p,ti)', ...
    Aq,T,S, xf), ...
    t, x0); 

subplot(3,3,7);
scatter(xout(:,1),xout(:,2),[],c); set(gca,'xlim',[-1 3]*1.2,'ylim',[-1 3]*1.2), colormap(gca,'jet')
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Line attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[8 9]);
plot(t,[v xout]);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')


function yd = odeAttractor2D( y, x, A, T, Gleak, fp)

    y = [y;1]; % make it homogeneous
    
    leak = y(1:2)'*[0 -1;1 0]*fp'; % TODO: need to generalize this cross product and i think and do Ay and Afp instead of the raw points.
    
    yd = - (4*y'*A*y) * T(:,:,1) * A*y ... % movement towards the attractor
                  + x * T(:,:,2) * A*y ... % movement along the attractor driven by the input
         - Gleak*leak * T(:,:,2) * A*y;... % movement along the attractor towards the fixed point
    
    yd = yd(1:end-1);

end

