%% attractor simulation
clear params;
close all


% time parameters of the simulation
tstart = 0;
tend = 10;
dt = 0.001;
trange = (tstart:dt:tend)';

% initial conditions of the attractor
e0 = [1-cosd(45) 1-sind(45)]';

% input velocity
x = zeros(size(trange));
x(trange > 1 & trange <2) = 1;
x(trange > 3 & trange <4) = 2;
x(trange > 5 & trange <9) = -1;
x = x+ randn(size(x))*.1;


% ring attractor
params.Wi = [0 -1 1; 1 0 -1]; % the first two columns act as a rotation matrix, the last column controls the center (should be minus the sum of the other numbers)
params.Wr = [0 0 0 ; 0 0 0]; % recurrent matrix BIG TODO. need to add some 
                            % dynamics to make the ring also able to drift and 
                            % to actually attract activity to the ring when it is not there
                            % but it is also ok to leave that problem to
                            % the decoder as the integrator will indeed
                            % integrate in rings
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
Orth1 = [1 1;1 -1]*(sqrt(2)/2); % decomposition of the weight matrix
D = diag([-1  -0.01]); % eigenvalues of the matrix

params.Wi = [0 0 .1; 0 0 -.1];
W = Orth1'*D*Orth1;
params.Wr = [W -sum(W,2)+e0-1]; % add a column at the end to shift the center of the line

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
    yd = params.Wi*yh*x  + params.Wr*yh;
end


