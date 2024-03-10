%% attractor simulation
clear params;
close all

% time parameters of the simulation
tstart = 0;
tend = 20;
dt = 0.001;
trange = (tstart:dt:tend)';

% initial conditions of the attractor
y0 = [1 0 0 0 0 0 0 0 0 0 0 0]';
y0 = [1 0 0 0 ]';
 y0 = [1 0 0 ]';

% input velocity
x = zeros(size(trange));
% x(trange > 1 & trange <2) = 10;
% x(trange > 3 & trange <4) = 20;
x(trange > 5 & trange <9) = -10;
x = x+ randn(size(x))*.1;

% ring attractor

A = 1;
B = 0;
C = 1;
D = 0;
E = 0;
F = -1;
params.Aq  = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];


n = length(y0);
phases =  pi + ((1:n)-1) * 2*pi ./ n;
p1 = (1+cos(0 + phases)+1+cos(pi + phases))/2;
p2 = (1+cos(0 + phases));
p3 = (1+cos(pi/2 + phases));
a1 = p2-p1;
a2 = p3-p1;

%     b2 = a2-(a2*a1')/(a1*a1')*a1;
%     way to get a coplanar vector to a1
%     and a2 that is orthogonal to a1;

B = zeros(3,n+1);
B(1,1:end-1) = a1/norm(a1);
B(2,1:end-1) = a2/norm(a1);
B(end,1:end-1) = p1;
B(end,end)=1;


params.FP = zeros(n,1); % no fixed point for now.
params.Gleak = 0*0.2;
params.B = B;

[t, e] = ode45(@(t,y)ode2DAttractor(t,y,trange,x,params), trange, y0);

figure
subplot(3,3,1);
plot(e(:,1),e(:,2)); set(gca,'xlim',[-1 3],'ylim',[-1 3])
xlabel( 'Internal unit 1'), ylabel( 'Internal unit 2')
title('Ring attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[2 3]);
% plot(t,[x exp(length(y0)*e)]);
plot(t,[x e]);
legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')
subplot(3,3,[4:9]);
imagesc(e')




figure

C= cov(e);
M= mean(e);
[V,D]= eig(C);
P = V * diag(sqrt(1./(diag(D) + 0.1))) ;
W1 = bsxfun(@minus, e, M);
W = W1 * P;

subplot(1,length(y0), 1)
plot(W(:,end),W(:,end-1),'o')

for i=1:length(y0)-1
    subplot(1,length(y0), i+1)
    plot(e(:,i),e(:,i+1),'o')
end

function yd = ode2DAttractor(t,y,xt,x,params)
    x = interp1(xt,x,t); % Interpolate the data set (gt,g) at time t
    

    y = [y;1]; % make it homogeneous

    % x is the input
    % y is the state/output of the atractor

    W = params.Aq;
    Ortho1 = [0 -1 0;1 0 0; 0 0 0]; %rotates Wy to get the orthogonal direction along the attractor
    Gleak = params.Gleak; % leak speed

    fp = params.FP;
    leak = y(2)*fp(1)-y(1)*fp(2);
    
    y = params.B*y;
    yd = (x-Gleak*leak)*Ortho1*W*y   -   (4*y'*W*y)*W*y;

    yd = params.B'*yd;
    yd = yd(1:end-1);
    
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

