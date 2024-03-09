% CONIC POLAR FORM
close all

theta = deg2rad(0:0.1:360);


es = [0.01 0.1 0.7 1 2 10000];

figure
hold
rotAngle = deg2rad(45);
for i=1:length(es)
    ep = es(i);
    r = ep ./(1-ep*cos(theta-rotAngle));

    plot(r.*cos(theta) +ep ./(1+ep)*cos(rotAngle) +1, r.*sin(theta)+ep ./(1+ep)*sin(rotAngle) +1,'o-')
end

set(gca,'xlim',[-0 3], 'ylim', [-0 3])

%%
close all

[x1,x2] = meshgrid(-1:0.1:2,-1:0.1:2);

X = [x1(:) x2(:) ones(size(x1(:)))];

%circle

A = 1;
B = 0;
C = 1;
D = -2;
E = -2;
F = 1;

% % line
% 
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
% 
% A = 1;
% B = 0;
% C = -1;
% D = 0;
% E = 0;
% F = -1;



Aq = [A B/2 D/2; B/2 C E/2; E/2 D/2 F ];


y = zeros(size(x1(:)));
for i=1:height(X)
    y(i) = X(i,:)*Aq*X(i,:)'*X(i,:)*Aq*X(i,:)';
end


g = zeros(height(x1(:)),3);
for i=1:height(X)
    g(i,:) = 4*X(i,:)*Aq*X(i,:)'*Aq*X(i,:)';
end

figure
contour(x1,x2,reshape(y,size(x1)),40)
hold
quiver(x1,x2,-reshape(g(:,1),size(x1)), -reshape(g(:,2),size(x1)))