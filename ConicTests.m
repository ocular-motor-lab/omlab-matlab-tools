% CONIC POLAR FORM
close all

theta = deg2rad(0:0.1:360);


es = [ 0.1 0.4 0.6 0.8 1 2 5 100];

figure
hold
rotAngle = deg2rad(45);
for i=1:length(es)
    ep = es(i);
    r = 1 ./(1-ep*cos(theta-rotAngle));
r(ep*cos(theta-rotAngle)>1) = nan;
    plot(r.*cos(theta) +1 ./(1+ep)*cos(rotAngle) +4, r.*sin(theta)+1 ./(1+ep)*sin(rotAngle) +4,'-')
end

set(gca,'xlim',[-0 10], 'ylim', [-0 10])

%%
close all

[x1,x2] = meshgrid(-1:0.1:3,-1:0.1:3);

X = [x1(:) x2(:) ones(size(x1(:)))];

%circle

% A = 1;
% B = 0;
% C = 1;
% D = -2;
% E = -2;
% F = 1;

% % line
% % 
A = 0;
B = 0;
C = 0;
D = 1;
E = 1;
F = -2;

% 
% % Ellipse
% 
% A = 1;
% B = 0;
% C = 5;
% D = 0;
% E = 0;
% F = -1;

% 
% 
% % hyperbola

% A = 1;
% B = 0;
% C = -1;
% D = 0;
% E = 0;
% F = -1;

% % parabola

A = 0;
B = -1;
C = 1;
D = -1;
E = -1;
F = 1;

Aq = [A B/2 D/2; B/2 C E/2; D/2 E/2 F ]     ;

y = zeros(size(x1(:)));
for i=1:height(X)
    y(i) = X(i,:)*Aq*X(i,:)'*X(i,:)*Aq*X(i,:)';
end


g = zeros(height(x1(:)),3);
for i=1:height(X)
    g(i,:) = 4*X(i,:)*Aq*X(i,:)'*Aq*X(i,:)';
end

figure
subplot(1,2,1)
contour(x1,x2,reshape(log(.1+y),size(x1)),40)
hold
quiver(x1,x2,-reshape(g(:,1),size(x1)), -reshape(g(:,2),size(x1)))
set(gca,'PlotBoxAspectRatio',[1 1 1])

d = zeros(height(x1(:)),3);
c = zeros(height(x1(:)),3);
for i=1:height(X)
    d(i,:) = [0 -1 0;1 0 0; 0 0 0]*Aq*X(i,:)';
    c(i,:) = Aq*X(i,:)';
end

subplot(1,2,2)
contour(x1,x2,reshape(log(.1+y),size(x1)),40)
hold
 quiver(x1,x2,-reshape(d(:,1),size(x1)), -reshape(d(:,2),size(x1)))
%quiver(x1,x2,-reshape(g(:,1),size(x1)), -reshape(g(:,2),size(x1)),'r')

set(gca,'PlotBoxAspectRatio',[1 1 1])

figure
% plot3(c(:,1),c(:,2),c(:,3),'o')
quiver(x1,x2,-reshape(c(:,1),size(x1)), -reshape(c(:,2),size(x1)),'r')