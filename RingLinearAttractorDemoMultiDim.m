%% attractor simulation
clear params;
close all

% time parameters of the simulation
dt = 0.001;
t = (0:dt:20)';

% initial conditions of the attractor
n = 20;
y0 = zeros(n,1);
y0(5) = 5;
% y0 = randn(n,1);
% input velocity
x = zeros(size(t));
x(t > 1 & t <2) = 10;
x(t > 3 & t <4) = 20;
x(t > 5 & t <9) = -15;
x = x+ randn(size(x))*0;
x=x*0.1;
% ring attractor
A = 1;
B = 0;
C = 1;
D = 0;
E = 0;
F = -2;
Aq  = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
T = cat(3, ...
    [1 0 0 ;    0 1 0 ;  0 0 0;], ...
    [0 -1 0 ;   1 0 0 ;  0 0 0;]);
% convoluted way to find two orthogonal vectors to be a base of the
% subspace that contains the manifold.
phases =  ((1:n)-1) * 2*pi ./ n;
p1 = zeros(size(phases));
p2 = cos(0 + phases);
p3 = cos(pi/2 + phases);
a1 = p2-p1;
a2 = p3-p1;
a1=a1/norm(a1);
a2=a2/norm(a2);
%     b2 = a2-(a2*a1')/(a1*a1')*a1;
%     way to get a coplanar vector to a1
%     and a2 that is orthogonal to a1;
S = [a1' a2']; % basis of the subspace plane containing the manifold
P = inv(S'*S)*S'; % projection matrix to the plane
P = [P zeros(2,1); zeros(1,n) 1]; % homogenous projection matrix to allow translation of the plane away from zero 
S = [S ones(n,1); zeros(1,2) 1]; 
[t, yout] = ode45(@(ti,yi)odeAttractorND(yi,interp1(t,x,ti), Aq, S, P, T), t, y0);
PlotRun(t,x,yout);
function yd = odeAttractorND( y, x, A, S, P, T)
    y = [y;1]; % make it homogeneous
    
    yb = P*y; % project into the plane that contains the manifold
    
    yd = - (4*yb'*A*yb) * T(:,:,1) * A*yb ... % movement towards the attractor
                  + x * T(:,:,2) * A*yb ;   % movement along the attractor driven by the input
    yd = S*yd + 10*(S*yb-y); % need to add this term to attract to the plane. TODO: not sure
    yd = yd(1:end-1);
end
function PlotRun(t,x,yout)
    figure
    subplot(3,3,1);
    c = linspace(1,255,length(yout));
    scatter(yout(:,1),yout(:,2),[],c); set(gca,'xlim',[-1 3]*1.2,'ylim',[-1 3]*1.2), colormap(gca,'jet')
    xlabel( 'unit 1'), ylabel( 'unit 2')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(3,3,4);
    scatter(yout(:,2),yout(:,3),[],c);  set(gca,'xlim',[-1 3]*1.2,'ylim',[-1 3]*1.2), colormap(gca,'jet')
    xlabel( 'unit 2'), ylabel( 'unit 3')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(3,3,7);
    scatter(yout(:,3),yout(:,4),[],c);  set(gca,'xlim',[-1 3]*1.2,'ylim',[-1 3]*1.2), colormap(gca,'jet')
    xlabel( 'unit 3'), ylabel( 'unit 4')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(3,3,[2 3]);
    plot(t,x,'linewidth',2);hold
    plot(t,yout);
    legend({'Input velocity x'})
    xlabel('Time')
    title('Wide tunning - cos(angle)')
    subplot(3,3,[5 6]);
    plot(t,x,'linewidth',2);hold
    plot(t,exp(width(yout)*yout));
    legend({'Input velocity x'})
    xlabel('Time')
    set(gca,'ylim',  exp(width(yout)*[min(min(yout(t>1,:))), max(max(yout(t>1,:)))]))
    title('Narrow tunning  - e^\(Nunits*cos(angle))')
    subplot(3,3,[8 9]);
    imagesc(exp(width(yout)*yout)');
    set(gca,'clim', exp(width(yout)*[min(min(yout(t>1,:))), max(max(yout(t>1,:)))])) % make the clim ignore the first 100 timepoints
    xlabel('Time')
    ylabel('units')
    set(gca,'xticklabel',[],'yticklabel',[])
end