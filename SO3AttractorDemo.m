%% SO3 attractor simulation
%% Jorge Otero-Millan 3/10/2024
% close all
N = 4; % 2 ring (SO1), 4 quaternions (SO3), 3 sphere O3 (not sure it will work) ;
n = 4;
rng(1);

% time parameters of the simulation
dt = 0.001;
t = (0:dt:20)';

% initial conditions of the attractor
% x0 = 2*[1 0 0 0 1 0 0 0  1 0 0 0  1 0 0 0  1 0 0 0]'; % point outside the attractor to show initial drift towards it
x0 = randn(n,1); % point outside the attractor to show initial drift towards it
x0 = [1 0 0 0];

% velocity input
w = zeros(length(t),N-3);
w(t > 3 & t <5, 1) = deg2rad(40); % angular velocity around x
w(t > 7 & t <10, 2) = deg2rad(60); % angular velocity around y
w(t > 12 & t <16, 3) = deg2rad(80); % angular velocity around z

w = w(:,1:N-1);

% so3 attractor matrix definition, conic section extended to be the
% hypersphere containing SO3. x'Ax =0 if x1^2 + x2^2 + x3^2 + x4^2 - 1 = 0
A = eye(N+1);
A(end,end)=-1;

% base of the subspace that contains the attractor
if ( N==n)
    S = eye(N);
else
    S = sin(2*pi*repmat((0:n-1)',1,N)/n.*(1+repmat(1:N,n,1)/300));
    % S = rand(n,N);
    S = orth(S);
end

% A tensor represents the group action inducing moving towards the
% attractor (first row) and along the atractor (next rows).
% The bottom part is the commutator
% TODO: Clean up and leave closer to Lie
switch(N)
    case 2
        T = cat(3, ...
            [0 -1 0 ;   1 0 0 ;  0 0 0;]);
    case 3
        T = cat(3, ...
            [0 -1 0 0; 0 0 1 0; -1 0 0 0; 0 0 0 0], ...
            [0 0 -1 0; 1 0 0 0; 0 1 0 0; 0 0 0 0]);
    case 4
        T = cat(3, ...
            [0 -1 0 0 0;   1 0 0 0 0;  0 0 0 -1 0; 0 0 1 0 0;  0 0 0 0 0], ...
            [0 0 -1 0 0;   0 0 0 1 0;  1 0 0 0 0;  0 -1 0 0 0; 0 0 0 0 0], ...
            [0 0 0 -1 0;   0 0 -1 0 0; 0 1 0 0 0;  1 0 0 0 0;  0 0 0 0 0]);
end

[t, xout] = ode45(@(ti,xi)odeAttractorND(...
    ti, xi, ...
    interp1(t,w,ti)', ...
    A,T,S), ...
    t, x0);

PlotRun(t,w,p, xout,inv(S'*S)*S'*xout');


%% 



function PlotRun(t,v,p,xout, xoutM)
    figure
    subplot(6,3,[1 4]);
    c = min(linspace(1,400,length(xout)),255);
    scatter(xout(:,1),xout(:,2),[],c); set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
    xlabel( 'unit 1'), ylabel( 'unit 2')
    title('SO3 attractor')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(6,3,[1 4]+6);
    scatter(xout(:,2),xout(:,3),[],c);  set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
    xlabel( 'unit 2'), ylabel( 'unit 3')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(6,3,[1 4]+12);
    scatter(xout(:,3),xout(:,4),[],c);  set(gca,'xlim',[-1 1]*1.2,'ylim',[-1 1]*1.2), colormap(gca,'jet')
    xlabel( 'unit 3'), ylabel( 'unit 4')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    subplot(6,3,[2 3]);
    plot(t,v,'linewidth',2);
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
    angles = rad2deg(quat2eul(xoutM','XYZ'));
    plot(t,angles);
    xlabel('Time')
    ylabel('angle(deg)')
    title('Decoded angle (euler XYZ)')

    h = get(gcf,'children');
    linkaxes(h([1 2 3 4 5]),'x');
end
