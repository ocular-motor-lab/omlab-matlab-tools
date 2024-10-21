

% %%
% clear all
% 
% tstart = 0;
% tend = 1;
% dt = 0.001;
% 
% % x(t) is a function that gets the value of x depending on t
% % TODO: it may be faster to consider t an index and prebuilt the x vector
% % then use dt to calculate time when needed. But this is more clear = @(t) ((t > 0.3 & t <0.4)*1 );
% trange = (tstart:dt:tend)';
% x = zeros(size(trange));
% x(trange > 0.3 & trange <0.4) = 1;
% 
% Tc = 0.1;
% params.W = 1/Tc;
% e0 = [0 1];
% 
% tic
% [t, e] = ode45(@(t,y)diffEq(t,y,trange,x,params), trange, e0);
% toc
% 
% tic
% % [t1, e1] = odeSimple(@(t,y)diffEq(t,y,trange,x,params), trange, e0);
% toc
% 
% % same but with filters
% Hc=tf([1],[Tc,1]);
% % Hc =
% %  
% %       1
% %   ---------
% %   0.1 s + 1
% %  
% % Continuous-time transfer function.
% Hd=c2d(Hc,dt);
% % Hd =
% %  
% %   0.00995
% %   --------
% %   z - 0.99
% %  
% % Sample time: 0.001 seconds
% % Discrete-time transfer function.
% b = Hd.Numerator{1};
% a = Hd.Denominator{1};
% 
% tic
% yf = filter(b,a,x);
% toc
% 
% 
% figure
% plot(t,[x,e],'-o');


%% Linear attractor simulation
clear params;

tstart = 0;
tend = 10;
dt = 0.001;
e0 = [0 1];

trange = (tstart:dt:tend)';
x = zeros(size(trange));
x(trange > 1 & trange <5) = 1;
x = sgolayfilt(x,1,11);

params.tau = 20;
params.V = [0 1]'*100;


% ring attractor
params.W0 =  zeros(2);
params.W1 =  [-1 0; 0 -1];
params.W2 =  [0 1; -1 0];
params.W3 =  [1 0; 0 1];
tic
[t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x,params), trange, e0);
toc

figure
subplot(1,2,1);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Ring attractor')
subplot(1,2,2,'nextplot','add');
plot(t,x);
plot(t,e);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')


% line attractor
Orth1 = [1 -1;1 1]*(sqrt(2)/2);
D = [-1 0;0 -0.01];

params.W0 =  Orth1'*D*Orth1;
params.W1 =  zeros(2);
params.W2 =  zeros(2);
params.W3 =  zeros(2);

tic
[t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x,params), trange, e0);
toc

figure
subplot(1,2,1);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Linear attractor')
subplot(1,2,2,'nextplot','add');
plot(t,x);
plot(t,e);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')

%%

% f(t,y)
function yd = diffEq(t,y,xt,x,params)
    x = interp1(xt,x,t); % Interpolate the data set (gt,g) at time t
    W = params.W;
%     yd = -y*W + x; integrator
    yd = -y*W + x*W;
end


function [t, yout] = odeSimple(fun, trange, y0)
if ( length(trange) == 2 )
    trange = linspace(trange(1),trange(2),1000);
end
yout= zeros(size(trange));
t = trange;

yout(1) = y0;
for i=2:length(t)
    yout(i) = yout(i-1) + fun(t(i), yout(i-1))*(t(i)-t(i-1)) ;
end

end


function yd = ode2DAttractor(t,y,xt,x,params)
    x = interp1(xt,x,t); % Interpolate the data set (gt,g) at time t

    if ( x>0)
        a=1;
    end
    v = x*params.V;
    y = y-[ 1 1]'; % shift center;

    yd = (params.W0)*v + (params.W1)*y(1)*v + (params.W2)*y(2)*v + (params.W3)*(y/(y'*y)-y);
    yd = yd/params.tau;
end