%% attractor simulation
clear params;
close all

tstart = 0;
tend = 10;
dt = 0.001;
e0 = [0 1];

trange = (tstart:dt:tend)';
x = zeros(size(trange));
x(trange > 1 & trange <5) = 1;
x = sgolayfilt(x,1,11) + randn(size(x))/10;

params.tau = 2;
params.V = [0 1]'*120;


% ring attractor
params.W0 =  zeros(2); % Matrix that pushes the state in a fixed direction depending on the input
params.W1 =  [-1 0; 0 -1]; % W1 and W2 behave as a rotating matrix between input an state
params.W2 =  [0 1; -1 0]; 
params.W3 =  [1 0; 0 1]; % Matrix to push the state towards the ring attractor of radious 1
params.W4 =  zeros(2);
[t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x,params), trange, e0);

figure
subplot(2,2,1);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Ring attractor')
subplot(2,2,2,'nextplot','add');
plot(t,x);
plot(t,e);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')


% line attractor
Orth1 = [1 -1;1 1]*(sqrt(2)/2);
D = diag([-1  -0.01]);

params.W0 =  eye(2);
params.W1 =  zeros(2);
params.W2 =  zeros(2);
params.W3 =  zeros(2);
params.W4 = Orth1'*D*Orth1;

[t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x,params), trange, e0);

subplot(2,2,3);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Linear attractor')
subplot(2,2,4,'nextplot','add');
plot(t,x);
plot(t,e);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')

%%
function yd = ode2DAttractor(t,y,xt,x,params)
    x = interp1(xt,x,t); % Interpolate the data set (gt,g) at time t
    
    v = x*params.V;
    y = y-[ 1 1]'; % shift center;

    % v is the input
    % y is the current state shifted
    yd = (params.W0)*v*(1) + (params.W1)*y(1)*(v) + (params.W2)*y(2)*(v)+ params.W4*y + (params.W3)*(y/(y'*y)-y) ;
    yd = yd/params.tau;
end