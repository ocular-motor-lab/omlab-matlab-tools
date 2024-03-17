%%
%% Basic conic shapes
close all

theta = deg2rad(0:0.1:360);


es = [ 0.1 0.4 0.6 0.8 1 2 5 100];

figure('color','w')
hold
rotAngle = deg2rad(45);
for i=1:length(es)
    ep = es(i);
    r = 1 ./(1-ep*cos(theta-rotAngle));
r(ep*cos(theta-rotAngle)>1) = nan;
    plot(r.*cos(theta) +1 ./(1+ep)*cos(rotAngle) +3, r.*sin(theta)+1 ./(1+ep)*sin(rotAngle) +3,'-','linewidth',2)
end
set(gca,'PlotBoxAspectRatio',[1 1 1],'fontsize',14)
set(gca,'xlim',[-0 10], 'ylim', [-0 10])
xlabel('x')
ylabel('y')

%% SO3 Attractor gradients

close all
theta = deg2rad(0:1:360);

% Circle quadratic
m = 2;
A = eye(3); A(end,end)=-1;

lim = 3;
[x1,x2] = meshgrid(-lim:0.1:lim, -lim:0.1:lim);
p = length(x2(:));
x = [x1(:),x2(:),ones(p,1)];

Ax = zeros(p,3);
xAx = zeros(p,1);
xAxAx = zeros(p,3);
for i=1:p
    Ax(i,:) = A*x(i,:)';
    xAx(i) = x(i,:)*A*x(i,:)';
    xAxAx(i,:) = 4*x(i,:)*A*x(i,:)'*A*x(i,:)';
end
xAxxAx = xAx.*(xAx);

figure('color','w')

subplot(2,2,1)
surf(x1, x2, reshape(xAx,size(x1)), 'facecolor',[1 1 1]*0.8,'EdgeColor', [1 1 1]*0.4)
zlim([-2 4])
ylim([-1 3])
xlim([-2 2])
hold
plot3(x1(:), x2(:), zeros(size(x1(:))),'.');
plot3(cos(theta),sin(theta), zeros(size(theta)),'linewidth', 2)
title('$q(x)=x^TAx$','interpreter','latex','fontsize',20)
xlabel('$x_1$','interpreter','latex','fontsize',16)
ylabel('$x_2$','interpreter','latex','fontsize',16)
set(gca,'xticklabel', [],  'yticklabel', [] )
set(gca,'zticklabel',{'' '0' ''})

subplot(2,2,2)
surf(x1, x2, reshape(xAx.*(xAx),size(x1)), 'facecolor',[1 1 1]*0.8,'EdgeColor', [1 1 1]*0.4)
zlim([-2 4])
ylim([-1 3])
xlim([-2 2])
hold
plot3(x1(:), x2(:), zeros(size(x1(:))),'.');
plot3(cos(theta),sin(theta), zeros(size(theta)),'linewidth', 2)
title('$q(x)^2$','interpreter','latex','fontsize',20)
xlabel('$x_1$','interpreter','latex','fontsize',16)
ylabel('$x_2$','interpreter','latex','fontsize',16)
set(gca,'xticklabel', [],  'yticklabel', [] )
set(gca,'zticklabel',{'' '0' ''})

subplot(2,2,3)
% surf(x1, x2, reshape(xAx.*(xAx),size(x1)), 'facecolor',[1 1 1]*0.9,'FaceAlpha',0.5, 'EdgeColor', [1 1 1]*0.8)

g = reshape(xAxAx, [size(x1),3]);
y = reshape(xAxxAx, size(x1));
% contour(x1,x2,y)
hold
idx = 1:4:size(x1,1);
quiver(x1(idx,idx),x2(idx,idx), -g(idx,idx,1),...
    -g(idx,idx,2),'linewidth',1)
plot(cos(theta),sin(theta), 'linewidth', 2)
% set(gca,'PlotBoxAspectRatio',[1 1 1])
title('$\nabla q(x)^2$','interpreter','latex','fontsize',20)
% view(-45,25)
zlim([-2 4])
ylim([-3 3])
xlim([-3 3])
xlabel('$x_1$','interpreter','latex','fontsize',16)
ylabel('$x_2$','interpreter','latex','fontsize',16)
set(gca,'xticklabel', [],  'yticklabel', [] )
set(gca,'zticklabel',[])

% set(gca, 'Clipping', 'off');

subplot(2,2,4)
theta2 = deg2rad(0:10:360);
x1 = cos(theta2)';
x2 = sin(theta2)';
p = length(theta2);
x = [x1(:),x2(:),ones(p,1)];

Ax = zeros(p,3);
xAx = zeros(p,1);
xAxAx = zeros(p,3);
L = zeros(p,2);
for i=1:p
    Ax(i,:) = A*x(i,:)';
    xAx(i) = x(i,:)*A*x(i,:)';
    xAxAx(i,:) = 4*x(i,:)*A*x(i,:)'*A*x(i,:)';
    L(i,:) = Ax(i,:)*[0 1; -1 0; 0 0];
end
xAxxAx = xAx.*(xAx);

g = Ax;%reshape(Ax, [size(x1),3]);
y = xAx;%reshape(xAx, size(x1));
% contour(x1,x2,y)
hold
idx = 1:4:size(x1,1);
quiver(x1,x2,-L(:,1), -L(:,2),1.3)
% set(gca,'PlotBoxAspectRatio',[1 1 1])
plot(cos(theta),sin(theta), 'linewidth', 2)
quiver(x1,x2,-g(:,1), -g(:,2),'color',[1 1 1]*0.8)
xlim([-1 1]*1.5)
ylim([-1 1]*1.5)
title('$ T  \nabla q(x)$','interpreter','latex','fontsize',20)
xlabel('$x_1$','interpreter','latex','fontsize',16)
ylabel('$x_2$','interpreter','latex','fontsize',16)
set(gca,'PlotBoxAspectRatio',[1 1 1])
set(gca,'xticklabel', [],  'yticklabel', [] , 'zticklabel', [])

%%
%%
%% Demo comparing sparse and non sparse coding
% sparse tunnin is just exp(N*x) where x is the non sparse tunning and N
% the number of units

close all

angle = 0:0.01:(2*pi);

Nunits = [2 3 5 8 20]; % number of units
NN = length(Nunits);

Phases = cell(size(Nunits));           % phases of each unit

X = cell(size(Nunits));             % activation of all units    
Xs = cell(size(Nunits));       % sparse activation of all units

for inu=1:NN
    X{inu} = zeros(length(angle), Nunits(inu));
    Phases{inu} = zeros(length(Nunits(inu)),1);

    for j=1:Nunits(inu)
        if ( Nunits(inu)>2)
            % the phases of each unit are distributed evenly along the circle
            Phases{inu}(j) = pi + (j-1) * 2*pi / Nunits(inu);
        else
            % the 2 unit case is special
            % units need to be in cuadrature otherwise they are not
            % decodable
            Phases{inu}(j) = (j-1) * pi/2;
        end
        
        % the activation of each unit is just a sinusoid, shifted by the
        % phase and elevated so it is always positive. Not critical but
        % more realistic if talking about neurons
        X{inu}(:,j) = 1 + cos( angle + Phases{inu}(j) );
    end
   

    % the sparse coding is just the exponential of the normal coding
    % the exponential includes a multiplicative term that makes the tuning
    % narrower the more units are present
    Xs{inu} = exp( Nunits(inu) * X{inu} ) / exp(Nunits(inu)*2);
end


NN = length(Phases);
figure('color','w')
for i= 1:NN
    subplot(5,NN,i,'nextplot','add')
    title(sprintf('%d units', Nunits(i)))
    %title('units and mean of units')
    plot(X{i})
  %  plot(mean(X{i},2),'linewidth',2,'color','k');
    set(gca,'xticklabel',[],'yticklabel',[])
    xlabel('Angle');

    subplot(5,NN,i+1*NN,'nextplot','add')
    %title('units and mean of units')
    plot(Xs{i})
   % plot(mean(Xs{i},2),'linewidth',2,'color','k');
    set(gca,'xticklabel',[],'yticklabel',[])
    xlabel('Angle');

    subplot(5,NN,i+2*NN)
    imagesc(Xs{i}');
    set(gca,'xticklabel',[],'yticklabel',[])
    ylabel('units');
    xlabel('Angle')

%     subplot(6,NN,i+3*NN)
%     plot(-atan2d(mean((X{i}-1).*repmat(sin(Phases{i}),height(X{i}),1),2), mean((X{i}-1).*repmat(cos(Phases{i}),height(X{i}),1),2)))  
%     hold
%     plot(-atan2d(mean((Xs{i}-1).*repmat(sin(Phases{i}),height(Xs{i}),1),2), mean((Xs{i}-1).*repmat(cos(Phases{i}),height(Xs{i}),1),2)))  
%     title('decoded angle')


    C= cov(X{i});
    M= mean(X{i});
    [V,D]= eig(C);
    P = V * diag(sqrt(1./(diag(D) + 0.1))) ;
    W1 = bsxfun(@minus, X{i}, M);
    W = W1 * P;


    C= cov(Xs{i});
    M= mean(Xs{i});
    [V,D]= eig(C);
    P = V * diag(sqrt(1./(diag(D) + 0.1))) ;
    W1 = bsxfun(@minus, Xs{i}, M);
    Ws = W1 * P;

    subplot(5,NN,i+3*NN)
        plot(W(:,end),W(:,end-1),'o')
        hold
        plot(Ws(:,end),Ws(:,end-1),'o')
        set(gca,'xlim',[-2 2],'ylim',[-2 2])
        set(gca,'xticklabel',[],'yticklabel',[])
        set(gca,'PlotBoxAspectRatio',[1 1 1])
        xlabel('PCA1');
        ylabel('PCA2');

    if (Nunits(i)>2)
    subplot(5,NN,i+4*NN)
        plot(W(:,end-1),W(:,end-2),'o')
        hold
        plot(Ws(:,end-1),Ws(:,end-2),'o')
        set(gca,'xlim',[-2 2],'ylim',[-2 2])
        set(gca,'xticklabel',[],'yticklabel',[])
        set(gca,'PlotBoxAspectRatio',[1 1 1])
        xlabel('PCA2');
        ylabel('PCA3');
    end
end

%%
% Change in shape of 2D manifold embedded in 3D with different tunings
Tunnings = logspace(1.5,-0.5,10);
X = {};
for k=1:length(Tunnings)

    Tunning=Tunnings(k);

    angle = 0:0.01:(2*pi);
    X{k} = zeros(length(angle), 3);
    for j=1:3
        Phase = (j-1)*2*pi/3;
        X{k}(:,j)= exp(sin(angle+Phase)/Tunning);
        X{k}(:,j) = X{k}(:,j)/mean(X{k}(:,j));
    end
end

figure ('color','white')
hold
for k=1:length(Tunnings)
    plot3(X{k}(:,1),X{k}(:,2),X{k}(:,3),'o')
end
view(135,45)


%% attractor simulation
clear params;
close all

% time parameters of the simulation
dt = 0.001;
t = (0:dt:20)';
c = t/max(t)*255;

% % line attractor
A = 0; B = 0; C = 0; D = 1; E = 1; F = -(D+E);
Aq  =  [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
T = cat(3, ...
    [0 -1;   1 0 ; ]);
S = eye(2);

% input velocity
w = zeros(size(t));
w(t > 1 & t <5) = deg2rad(90);
w(t > 8 & t <12) = deg2rad(-180);

% initial conditions of the attractor
x0 = [-1 -1]';

[t, xout] = ode45(@(ti,xi)AttractorNetwork(ti, xi, interp1(t,w/2,ti)', Aq,T,S), t, x0); 

figure('color','w')
subplot(3,3,1);
scatter(xout(:,1),xout(:,2),[],c); set(gca,'xlim',[-1 3]*1.2,'ylim',[-1 3]*1.2), colormap(gca,'jet')
xlabel( 'x_1'), ylabel( 'x_2')
title('Line attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[2 3]);
plot(t,[w xout]);
xlabel('Time')
legend({'\omega', 'x_1', 'x_2'})


% ring attractor
A = 1; B = 0; C = 1; D = 0; E = 0; F = -1;
Aq  = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
T = cat(3, ...
    [0 -1;   1 0 ; ]);
S = eye(2);

w = zeros(size(t));
w(t > 1 & t <2) = deg2rad(90);
w(t > 3 & t <4) = deg2rad(-90);
w(t > 8 & t <12) = deg2rad(180);
%v = v+ randn(size(v))*.1;
% initial conditions of the attractor
x0 = [-1 -1]';

[t, xout] = ode45(@(ti,xi)AttractorNetwork( ti, xi, interp1(t,w,ti)', Aq,T,S), t, x0); 

subplot(3,3,4);
scatter(xout(:,1),xout(:,2),[],c); set(gca,'xlim',[-1 1]*1.5,'ylim',[-1 1]*1.5), colormap(gca,'jet')
xlabel( 'x_1'), ylabel( 'x_2')
title('Ring attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[5 6]);
plot(t,[w xout]);
xlabel('Time')


% parabollic attractor
A = 0; B = -1; C = 1; D = -1; E = -1; F = 1;
Aq = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
T = cat(3, ...
    [0 -1;   1 0 ; ]);
S = eye(2);

w = zeros(size(t));
w(t > 1 & t <5) = deg2rad(90);
w(t > 8 & t <12) = deg2rad(-180);

[t, xout] = ode45(@(ti,xi)AttractorNetwork(  ti, xi,   interp1(t,w/4,ti)', Aq,T,S), t, x0); % velocity made slower to avoid numerical problems

subplot(3,3,7);
scatter(xout(:,1),xout(:,2),[],c); set(gca,'xlim',[-1 3]*1.2,'ylim',[-1 3]*1.2), colormap(gca,'jet')
xlabel( 'x_1'), ylabel( 'x_2')
title('Parabollic attractor')
set(gca,'PlotBoxAspectRatio',[1 1 1])
subplot(3,3,[8 9]);
plot(t,[w xout]);
% legend({'Input velocity', 'Internal unit 1', 'Internal unit 2'})
xlabel('Time')



%% SO3 attractor simulation
%% Jorge Otero-Millan 3/10/2024
% close all
N = 3; % 2 ring (SO1), 3 sphere (S2), 4 quaternions (SO3)
n = 6;
rng(1);

% time parameters of the simulation
dt = 0.001;
t = (0:dt:20)';


% velocity input
w = zeros(length(t),N-1);
rt = t(t > 3 & t <18);

% w(t > 3 & t <18, 1) = cos(rt*2)*1.5; % flower for sphere
% w(t > 3 & t <18, 2) = sin(rt*2)*1.5;

w(t > 3 & t <18, 1) = cos(rt*2)*1.5;
w(t > 3 & t <18, 2) = sin(rt*2)*1.5;

% w(t > 3 & t <18, 2) = deg2rad(180/1); % angular velocity around x
% w(t > 7 & t <18, 2) = deg2rad(180/4); % angular velocity around y
% w(t > 12 & t <16, 3) = deg2rad(75); % angular velocity around z
w = w(:,1:N-1); % get only the relevant componets

% so3 attractor matrix definition, conic section extended to be the
% hypersphere containing SO3. x'Ax =0 if x1^2 + x2^2 + x3^2 + x4^2 - 1 = 0
A = eye(N+1);
A(end,end)=-1;

% base of the subspace that contains the attractor
if ( N==n)
    S = eye(N);
else
    S = sin(2*pi*repmat((0:n-1)',1,N)/n.*(1+repmat(1:N,n,1)/300));
%     sin(mod(repmat((0:n-1)',1,N)/n+repmat((0:N-1),n,1)/4,1)*2*pi)/sqrt(2)
    % S = rand(n,N);
    S = orth(S);
end

S = S*sqrt(n-1)/sqrt(2);

% initial conditions of the attractor
x0 = S(:,1);

switch(N)
    case 2
        T = cat(3, ...
            [0 -1 ;  1 0 ]); % 2x2x1
    case 3
        T = cat(3, ...
            [0 -1  0 ;  1  0  0 ;  0  0  0 ], ...
            [0  0  1 ;  0  0  0 ; -1  0  0 ]); % 3x3x2
    case 4
        T = cat(3, ...
            [0 -1  0  0 ;  1  0  0  0 ;  0  0  0 -1 ; 0  0  1  0 ], ...
            [0  0 -1  0 ;  0  0  0  1 ;  1  0  0  0 ; 0 -1  0  0 ], ...
            [0  0  0 -1 ;  0  0 -1  0 ;  0  1  0  0 ; 1  0  0  0 ]); % 4x4x3
end



[t, xout] = ode45(@(ti,xi) ...
    ...
    AttractorNetwork( ti, xi, interp1(t,w,ti)', A, T, S), ...
    ...
    t, x0);

xoutM = inv(S'*S)*S'*xout';

figure('color','w')
subplot(6,3,[1 4]);
c = min(linspace(1,400,length(xout)),255);
scatter(xout(:,1),xout(:,2),[],c); set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
xlabel( 'x_1'), ylabel( 'x_2')
%     title('SO3 attractor')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(6,3,[1 4]+6);
    scatter(xout(:,2),xout(:,3),[],c);  set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
xlabel( 'x_2'), ylabel( 'x_3')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    if ( size(xout,2)>3)
        subplot(6,3,[1 4]+12);
        scatter(xout(:,3),xout(:,4),[],c);  set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
        xlabel( 'x_3'), ylabel( 'x_4')
        set(gca,'PlotBoxAspectRatio',[1 1 1])
    end
    
    subplot(6,3,[2 3]);
    plot(t,w,'linewidth',2);
    set(gca,'xticklabel',[],'yticklabel',[])
    title('Input velocity')

    subplot(6,3,[2 3 5 6]+3);
    plot(t,xout);
    hold
    plot(t, vecnorm(xout'),'k','linewidth',2)
%     set(gca,'xticklabel',[],'yticklabel',[])
    title('Network units')

    subplot(6,3,[2 3 5 6]+9);
    imagesc([t(1) t(end)],[1 width(xout)], xout')
    set(gca,'clim', [min(min(xout(t>1,:))), max(max(xout(t>1,:)))]) % make the clim ignore the first 100 timepoints
    ylabel('units')
    set(gca,'xticklabel',[],'yticklabel',[])

    subplot(6,3,[2 3]+15);
    angles = rad2deg(quat2eul(xoutM','XYZ'));
    plot(t,angles);
    xlabel('Time')
    ylabel('angle(deg)')
    title('Decoded angle (euler XYZ)')

    h = get(gcf,'children');
    linkaxes(h([1 2 3 4]),'x');

%% attractor simulation
clear params;
% close all

n = 20;

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
x0 = p1;
% y0 = randn(n,1);
% input velocity
w = zeros(size(t));
w(t > 1 & t <2) = 10;
w(t > 3 & t <4) = 20;
w(t > 5 & t <9) = -15;
w(t > 10 & t <20) = 45;
w = w+ randn(size(w))*0;
w=w*0.1;
% ring attractor
A = 1;
B = 0;
C = 1;
D = 0;
E = 0;
F = -1;
Aq  = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
T = cat(3, ...
    [0 -1 ;   1 0 ; ]);

% P = inv(S'*S)*S'; % projection matrix to the plane
% P = [P zeros(2,1); zeros(1,n) 1]; % homogenous projection matrix to allow translation of the plane away from zero 
% S = [S ones(n,1); zeros(1,2) 1]; 

p= zeros(length(t),2);

[t, xout] = ode45(@(ti,xi)AttractorNetwork(...
    ti, xi, ...
    interp1(t,w,ti)', ...
    Aq,T,S), ...
    t, x0); 

xoutM = inv(S'*S)*S'*xout';

    figure('color','w')
    subplot(6,3,[1 4]);
    c = min(linspace(1,400,length(xout)),255);
    scatter(xout(:,1),xout(:,2),[],c); set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
xlabel( 'x_1'), ylabel( 'x_2')
%     title('SO3 attractor')
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
    title('Decoded angle (deg)')

    h = get(gcf,'children');
    linkaxes(h([1 2 3 4 5]),'x');
