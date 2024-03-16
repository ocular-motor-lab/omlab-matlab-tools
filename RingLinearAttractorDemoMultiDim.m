%% attractor simulation
clear params;
% close all

% convoluted way to find two orthogonal vectors to be a base of the
% subspace that contains the manifold.
phases =  ((1:n)-1) * 2*pi ./ n;
p1 = cos(0 + phases);
p2 = cos(pi/2 + phases);
p1=p1/norm(p1);
p2=p2/norm(p2);
S = [p1' p2']; % basis of the subspace plane containing the manifold

% time parameters of the simulation
dt = 0.001;
t = (0:dt:20)';

% initial conditions of the attractor
n = 20;
x0 = p1;
% y0 = randn(n,1);
% input velocity
v = zeros(size(t));
v(t > 1 & t <2) = 10;
v(t > 3 & t <4) = 20;
v(t > 5 & t <9) = -15;
v(t > 10 & t <20) = 45;
v = v+ randn(size(v))*0;
v=v*0.1;
% ring attractor
A = 1;
B = 0;
C = 1;
D = 0;
E = 0;
F = -1;
Aq  = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
T = cat(3, ...
    [0 -1 0 ;   1 0 0 ;  0 0 0;]);

% P = inv(S'*S)*S'; % projection matrix to the plane
% P = [P zeros(2,1); zeros(1,n) 1]; % homogenous projection matrix to allow translation of the plane away from zero 
% S = [S ones(n,1); zeros(1,2) 1]; 

p= zeros(length(t),2);

[t, xout] = ode45(@(ti,xi)odeAttractorND(...
    ti, xi, ...
    interp1(t,v,ti)', ...
    Aq,T,S), ...
    t, x0); 

xoutM = inv(S'*S)*S'*xout';

    figure('color','w')
    subplot(6,3,[1 4]);
    c = min(linspace(1,400,length(xout)),255);
    scatter(xout(:,1),xout(:,2),[],c); set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
xlabel( 'x_1'), ylabel( 'x_2')
    title('SO3 attractor')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(6,3,[1 4]+6);
    scatter(xout(:,2),xout(:,3),[],c);  set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
xlabel( 'x_2'), ylabel( 'x_3')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(6,3,[1 4]+12);
    scatter(xout(:,3),xout(:,4),[],c);  set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
xlabel( 'x_3'), ylabel( 'x_4')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(6,3,[2 3]);
    plot(t,w,'linewidth',2);
    set(gca,'xticklabel',[],'yticklabel',[])
    title('Input velocity')

    subplot(6,3,[2 3]+3);
    plot(t,p);
    set(gca,'xticklabel',[],'yticklabel',[])
    title('Input position')

    subplot(6,3,[2 3]+6);
    plot(t,xout);
    set(gca,'xticklabel',[],'yticklabel',[])
    title('Network units')

    subplot(6,3,[2 3 5 6]+9);
    imagesc([t(1) t(end)],[1 width(xout)], xout')
    set(gca,'clim', [min(min(xout(t>1,:))), max(max(xout(t>1,:)))]) % make the clim ignore the first 100 timepoints
    ylabel('units')
    set(gca,'xticklabel',[],'yticklabel',[])

    subplot(6,3,[2 3]+15);
    angles =atan2d(xoutM(2,:),xoutM(1,:));
    plot(t,angles);
    xlabel('Time')
    ylabel('angle(deg)')
    title('Decoded angle (euler XYZ)')

    h = get(gcf,'children');
    linkaxes(h([1 2 3 4 5]),'x');