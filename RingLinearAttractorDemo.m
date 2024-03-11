%% attractor simulation
clear params;
close all

% time parameters of the simulation
dt = 0.001;
t = (0:dt:20)';

% initial conditions of the attractor
y0 = [0 0]';

% input velocity
x = zeros(size(t));
x(t > 1 & t <2) = 10;
x(t > 3 & t <4) = 20;
x(t > 5 & t <9) = -10;
x = x+ randn(size(x))*.1;

% ring attractor
A = 1;
B = 0;
C = 1;
D = -2;
E = -2;
F = 1;
Aq  = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
FP = [1 1] - [cosd(45) sind(45)];
T = cat(3, ...
    [1 0 0 ;    0 1 0 ;  0 0 0;], ...
    [0 -1 0 ;   1 0 0 ;  0 0 0;]);
Gleak = 0.2;

[t, yout] = ode45(@(ti,yi)odeAttractor2D(yi,interp1(t,x,ti),Aq,T,Gleak,FP), t, y0);

figure
    c = linspace(1,255,length(yout));
subplot(3,3,1);
scatter(yout(:,1),yout(:,2),[],c); set(gca,'xlim',[-1 3]*1.2,'ylim',[-1 3]*1.2), colormap(gca,'jet')
xlabel( 'Internal unit 1'), ylabel( 'Internal unit 2')
title('Ring attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[2 3]);
plot(t,[x yout]);
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
FP = -[0.5 0.5];

[t, yout] = ode45(@(ti,yi)odeAttractor2D(yi,interp1(t,x/10,ti),Aq,T,Gleak,FP), t, y0); % velocity made slower to avoid numerical problems

subplot(3,3,4);
scatter(yout(:,1),yout(:,2),[],c); set(gca,'xlim',[-1 3]*1.2,'ylim',[-1 3]*1.2), colormap(gca,'jet')
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Parabollic attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[5 6]);
plot(t,[x yout]);
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
FP = [1 1];

[t, yout] = ode45(@(ti,yi)odeAttractor2D(yi,interp1(t,x,ti),Aq,T,Gleak,FP), t, y0);

subplot(3,3,7);
scatter(yout(:,1),yout(:,2),[],c); set(gca,'xlim',[-1 3]*1.2,'ylim',[-1 3]*1.2), colormap(gca,'jet')
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Line attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[8 9]);
plot(t,[x yout]);
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

