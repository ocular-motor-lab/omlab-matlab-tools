%% attractor simulation
clear params;
close all

% time parameters of the simulation
tstart = 0;
tend = 10;
dt = 0.001;
trange = (tstart:dt:tend)';

% initial conditions of the attractor
e0 = [0 0]';

% input velocity
x = zeros(size(trange));
x(trange > 1 & trange <2) = 10;
x(trange > 3 & trange <4) = 20;
x(trange > 5 & trange <9) = -10;
x = x+ randn(size(x))*.1;

e0 = [0 0]';

% ring attractor
params.Wi = [0 -1 1; 1 0 -1]; 
% the first two columns act as a rotation matrix, the last column controls
% the center (should be minus the sum of the other numbers) 

A = 1;
B = 0;
C = 1;
D = -2;
E = -2;
F = 1;
params.Wr  = [A B/2 D/2; B/2 C E/2; E/2 D/2 F ];

% conic representation of the attractor

[t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x,params), trange, e0);

figure
subplot(2,2,1);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
xlabel( 'Internal unit 1'), ylabel( 'Internal unit 2')
title('Ring attractor')
subplot(2,2,2);
plot(t,[x e]);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')


% line attractor
A = 0;
B = 0;
C = 0;
D = 1;
E = 1;
F = -(D+E);
params.Wr  = [A B/2 D/2; B/2 C E/2; E/2 D/2 F ];

params.Wi = [0 0 .1; 0 0 -.1];

[t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x,params), trange, e0);

subplot(2,2,3);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Line attractor')
subplot(2,2,4);
plot(t,[x e]);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')


function yd = ode2DAttractor(t,y,xt,x,params)
    x = interp1(xt,x,t); % Interpolate the data set (gt,g) at time t
    
    yh = [y;1]; % make it homogeneous

    % x is the input
    % y is the state

    % the first term is the effect of the input non linearly associated
    % with the current state to allow for rotations

    % the second term is the one that attracts to the conic curve
    % which can be a line, a circle or whatever. 

    % right now missing is the additional point attractor to allow for
    % leakiness. 
    yd = params.Wi*yh*x - [1 0 0 ; 0 1 0]*(4*yh'*params.Wr*yh*params.Wr*yh);
end


