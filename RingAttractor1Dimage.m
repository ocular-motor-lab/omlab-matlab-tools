%% attractor simulation
clear params;
 % close all
% time parameters of the simulation
dt = 0.01;
t = (0:dt:20)';

% initial conditions of the attractor
n = 3;
k=200;

I = randn(k,1);

phases =  (0:n-1) * 2*pi ./ n;
p1 = cos(0 + phases);
p2 = sin(0 + phases);
S = [p1' p2'];
x0=p2;

phases =  (0:k-1) * pi ./ k;
p1 = cos(0 + phases);
p2 = sin(0 + phases);
R = [p1' p2'];

RS = R*S';


% input velocity
v = zeros(size(t));
% v(t > 0 & t <=1) = deg2rad(90);
v(t > 4 & t <=16) = -deg2rad(90);
% v(t > 5 & t <9) = -15;

% ring attractor
A = 1;
B = 0;
C = 1;
D = 0;
E = 0;
F = -1;
Aq  = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ];
T = cat(3, ...
    [1 0 0 ;    0 1 0 ;  0 0 0;], ...
    [0 -1 0 ;   1 0 0 ;  0 0 0;]);
% convoluted way to find two orthogonal vectors to be a base of the
% subspace that contains the manifold.
phases =  (0:n-1) * 2*pi ./ n;
p1 = zeros(size(phases));
p2 = cos(0 + phases);
p3 = sin(0 + phases);
% a1 = p2-p1;
% a2 = p3-p1;
% a1=a1/norm(a1);
% a2=a2/norm(a2);
%     b2 = a2-(a2*a1')/(a1*a1')*a1;
%     way to get a coplanar vector to a1
%     and a2 that is orthogonal to a1;
% S = [a1' a2']; % basis of the subspace plane containing the manifold
S = [p2' p3'];
% S=randn(size(S));
% P = inv(S'*S)*S'; % projection matrix to the plane
% P = [P zeros(2,1); zeros(1,n) 1]; % homogenous projection matrix to allow translation of the plane away from zero 
% S = [S ones(n,1); zeros(1,2) 1]; 

p= zeros(length(t),2);

[t, xout] = ode45(@(ti,xi)odeAttractorImage(...
    ti, xi, ...
    interp1(t,v,ti)', interp1(t,p,ti)', ...
    Aq,T,S), ...
    t, x0); 

[n, m] = size(S); % dimensions space and subspace
P = inv(S'*S)*S';
P = [P zeros(m,1); zeros(1,n) 1]; 
xm = (P*[xout ones(height(xout),1)]')';
S = [S zeros(n,1); zeros(1,m) 1];

PlotRun(t,v,xout, xm(:,1:end-1),S);


Iout = zeros(length(t),length(I));
for i=1:length(t)
   Iout(i,:) = RS*diag(xout(i,:))*RS'*I;
end


figure
imagesc(Iout)

function xd = odeAttractorImage( ti, xi, vi, pi, A, T, S)
% ODEATTRACTORND  differential equation for an attractor network of m
% dimensions embedded in an n dimensional space with p inputs. 
%   xd = odeAttractorND( ti, xi, vi, pi, A, T, S, xf) 
%
%   outputs:
%       xd: derivative of the state x at time ti
%
%   inputs: 
%       ti: (1 x 1) time 
%       xi: (n x 1) state at time ti 
%       vi: (q x 1) velocity input at time ti 
%       pi: (n x 1) position input at time ti 
%       A:  (m+1 x m+1) matrix defining the quadratic ecuation that specifis
%           the attractor shape. The attractor is the set of points x that
%           meet the condition x'*A*x=0
%       T: (m+1 x m+1 x q+1) Tensor encoding the directions of
%           motion from x relative to the direction orthogonal to the
%           attractor. Each of the p+1 matrices roates the direction
%           towards the attractor so we get one direction typically towards
%           the attractor and p directions along the attractor.
%       S: (n x m) Matrix to project x into a subspace that contains the
%           attractor. Each column is a base vector defining the subspace.
%       xf: (n x 1) Fixed point (set point or null point) on the space. The
%           state will drift towards it in the absence of inputs.
%

    [n, m] = size(S); % dimensions space and subspace
    q = size(T,3)-1; % number of input dimensions
    
    P = inv(S'*S)*S'; % projection matrix to the plane if S not
    % orthonormal
    P = [P zeros(m,1); zeros(1,n) 1]; % homogenous projection matrix to allow translation of the plane away from zero
    S = [S zeros(n,1); zeros(1,m) 1];

    x  = [xi;1]; % make it homogeneous
    p  = [pi;1];
    v  = vi;
    
    % Tensor product to rotate an orthogonal frame and align it with the
    % surface of the attractor.
    % Dim1 is the direction of movement towards the attractor
    % Dims 1..n are the directions of movement along the attractor dim1
    % driven by the input 
    M = squeeze(pagemtimes(T, A*P*x));

    % Scale the directions of motion by the velocity
    xd = M * [ - 0*(4*x'*P'*A*P*x);    ...    % velocity towards the attractor
              v ];                          % velocity along the attractor scaled by velocity input

    % remove the homogenous component
    xd = S*xd;
    xd = xd(1:end-1);
end

% IDEA for the drifts

% if  Ax and (xf-x) are parallel it means that we don't need to move

function PlotRun(t,x,yout, ym, S)
    figure
    subplot(3,3,1,'nextplot','add');
    c = linspace(1,255,length(yout));
    scatter(ym(:,1), ym(:,2),[],c); 
%     set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2); 
    colormap(gca,'jet')
    xlabel( 'pca 1'), ylabel( 'pca 2')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
        line([0 0],get(gca,"YLim"));
        line(get(gca,"XLim"),[0 0]);

    if ( width(yout)>1)
        subplot(3,3,4,'nextplot','add');
        c = linspace(1,255,length(yout));
        scatter(yout(:,1),yout(:,2),[],c);
%         set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2);
        colormap(gca,'jet')
        xlabel( 'unit 1'), ylabel( 'unit 2')
        set(gca,'PlotBoxAspectRatio',[1 1 1])
        line([0 0],get(gca,"YLim"));
        line(get(gca,"XLim"),[0 0]);
    end

    if ( width(yout)>2)
        subplot(3,3,7, 'nextplot','add');
        scatter3(yout(:,1),yout(:,2),yout(:,3));  
%         set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), 
        colormap(gca,'jet')
        xlabel( 'unit 2'), ylabel( 'unit 3')
        set(gca,'PlotBoxAspectRatio',[1 1 1])
        line([0 0],[-6 6], [0 0]);
        line([-6 6],[0 0], [0 0]);
        line([0 0], [0 0], [-6 6]);
%         line(get(gca,"XLim"),[0 0]);

        [X,Y] = meshgrid([-2:0.2:2],[-2:0.2:2]);
        pp = (S(1:end-1,1:end-1)*[X(:) Y(:)]')';
        c2 = linspace(1,255,length(pp));
        scatter3(pp(:,1),pp(:,2),pp(:,3),'.');  

    end
    
    subplot(3,3,[2 3]);
    plot(t,x,'linewidth',2);hold
    legend({'Input velocity x'})
    xlabel('Time')
    subplot(3,3,[5 6]);
    plot(t,yout);
    xlabel('Time')
%     set(gca,'ylim',  exp(width(yout)*[min(min(yout(t>1,:))), max(max(yout(t>1,:)))]))
%     title('Narrow tunning  - e^\(Nunits*cos(angle))')
    subplot(3,3,[8 9]);
    imagesc(yout');
%     imagesc(exp(width(yout)*yout)');
%     set(gca,'clim', exp(width(yout)*[min(min(yout(t>1,:))), max(max(yout(t>1,:)))])) % make the clim ignore the first 100 timepoints
    xlabel('Time')
    ylabel('units')
    set(gca,'xticklabel',[],'yticklabel',[])
end