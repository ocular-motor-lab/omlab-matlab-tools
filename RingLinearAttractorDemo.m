%% attractor simulation
clear params;
close all

tstart = 0;
tend = 10;
dt = 0.001;
e0 = [1-cosd(45) 1-sind(45)];

trange = (tstart:dt:tend)';
x = zeros(size(trange));
x(trange > 1 & trange <2) = 10;
x(trange > 3 & trange <4) = 10;
x(trange > 5 & trange <7) = -10;
x = x+ randn(size(x))*1;

params.tau = 1;


% ring attractor
params.W1 = [0 -1 1; 1 0 -1];
params.W2 = [0 0 0 ; 0 0 0];
[t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x,params), trange, e0);

figure
subplot(3,2,1);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Ring attractor')
subplot(3,2,2,'nextplot','add');
plot(t,x);
plot(t,e);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')


% line attractor
e0 = [1 1];
Orth1 = [1 1;1 -1]*(sqrt(2)/2);
D = diag([-1  -0.0001]);

params.W1 = [0 0 .1; 0 0 -.1];
Modes = Orth1'*D*Orth1;
params.W2 = [Modes -sum(Modes,2)]; % add a column at the end to shift the center of the line

[t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x,params), trange, e0);

subplot(3,2,3);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
xlabel( 'Internal unit 1')
ylabel( 'Internal unit 2')
title('Line attractor')
subplot(3,2,4,'nextplot','add');
plot(t,x);
plot(t,e);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')



% % line attractor
% e0 = [1 1];
% Orth1 = [1 1;1 -1]*(sqrt(2)/2);
% D = diag([-1  -0.1]);
% 
% params.W1 = 0.2*1/2*[0 -1 0; 1 0 0] + 1/2*[0 0 .1; 0 0 -.1];
% params.W2 = 0.2*1/2*[0 0; 0 0] + 1/2*Orth1'*D*Orth1;
% 
% [t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x,params), trange, e0);
% 
% subplot(3,2,5);
% plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
% xlabel( 'Internal unit 1')
% ylabel( 'Internal unit 2')
% title('Mixed attractor')
% subplot(3,2,6,'nextplot','add');
% plot(t,x);
% plot(t,e);
% legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
% xlabel('Time')



%%
% function yd = ode2DAttractor(t,y,xt,x,params)
%     x = interp1(xt,x,t); % Interpolate the data set (gt,g) at time t
%     
%     y = y-[ 1 1]'; % shift center;
% 
%     % v is the input
%     % y is the current state shifted
%     yd = (params.W0)*v*(1) + (params.W1)*y(1)*(v) + (params.W2)*y(2)*(v)+ params.W4*y + (params.W3)*(y/(y'*y)-y) ;
%     yd = yd/params.tau;
% end
%%
function yd = ode2DAttractor(t,y,xt,x,params)
    x = interp1(xt,x,t); % Interpolate the data set (gt,g) at time t
    
    yh = [y;1]; % make it homogeneous

    % x is the input
    % y is the state
    yd = params.W1*yh*x  + params.W2*yh; %  (params.W3)*(y/(y'*y)-y) ;
    yd = yd/params.tau;
end


% CONIC POLAR FORM
% close all
% 
% theta = deg2rad(0:0.1:360);
% 
% 
% es = [0.01 0.1 0.7 1 2 100];
% 
% figure
% hold
% for i=1:length(es)
%     ep = es(i);
%     r = ep ./(1-ep*cos(theta-pi/4));
% 
%     plot(r.*cos(theta)+1,r.*sin(theta)+1,'o-')
% end
% 
% set(gca,'xlim',[-0 3], 'ylim', [-0 3])