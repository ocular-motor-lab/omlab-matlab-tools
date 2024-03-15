%%
% SOn Attractor Figures

close all

% Circle quadratic
m = 2;
A = eye(3); A(end,end)=-1;

lim = 2;
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
surf(x1, x2, reshape(xAx,size(x1)));
zlim([-2 4])
ylim([-1 7])
hold
surf(x1, x2, zeros(size(x1)),'facecolor','red');
title('$x^TAx$','interpreter','latex','fontsize',16)

subplot(2,2,2)
surf(x1, x2, reshape(xAx.*(xAx),size(x1)));
zlim([-2 4])
ylim([-1 7])
hold
surf(x1, x2, zeros(size(x1)),'facecolor','red');
title('$x^TAxx^TAx$','interpreter','latex','fontsize',16)

subplot(2,2,3)
theta = deg2rad(0:1:360);

g = reshape(Ax, [size(x1),3]);
y = reshape(xAx, size(x1));
% contour(x1,x2,y)
hold
idx = 1:4:size(x1,1);
quiver(x1(idx,idx),x2(idx,idx),-g(idx,idx,1),...
    -g(idx,idx,2))
set(gca,'PlotBoxAspectRatio',[1 1 1])
title('$Ax$','interpreter','latex','fontsize',16)
plot(cos(theta),sin(theta), 'linewidth', 2)
subplot(2,2,4)

g = reshape(xAxAx, [size(x1),3]);
y = reshape(xAxxAx, size(x1));
% contour(x1,x2,y)
hold
idx = 1:2:size(x1,1);
quiver(x1(idx,idx),x2(idx,idx),-g(idx,idx,1),...
    -g(idx,idx,2))
set(gca,'PlotBoxAspectRatio',[1 1 1])
title('$4x^TAx^TAx$','interpreter','latex','fontsize',16)
plot(cos(theta),sin(theta), 'linewidth', 2)

%%
% Sphere quadratic
m = 2;
A = eye(4); A(end,end)=-1;

lim = 2;
[x1,x2, x3] = meshgrid(-lim:0.1:lim, -lim:0.1:lim,-lim:0.1:lim);
p = length(x2(:));
x = [x1(:),x2(:),x3(:),ones(p,1)];

Ax = zeros(p,4);
xAx = zeros(p,1);
xAxAx = zeros(p,4);
for i=1:p
    Ax(i,:) = A*x(i,:)';
    xAx(i) = x(i,:)*A*x(i,:)';
    xAxAx(i,:) = 4*x(i,:)*A*x(i,:)'*A*x(i,:)';
end
xAxxAx = xAx.*(xAx);

figure('color','w')
subplot(2,2,1)

% surf(reshape(x(1,:), size(x1)),reshape(x(1,:), size(x1),reshape(x(1,:), size(x1))

return


surf(x1, x2, reshape(xAx,size(x1)));
zlim([-2 4])
ylim([-1 7])
hold
surf(x1, x2, zeros(size(x1)),'facecolor','red');
title('$x^TAx$','interpreter','latex','fontsize',16)

subplot(2,2,2)
surf(x1, x2, reshape(xAx.*(xAx),size(x1)));
zlim([-2 4])
ylim([-1 7])
hold
surf(x1, x2, zeros(size(x1)),'facecolor','red');
title('$x^TAxx^TAx$','interpreter','latex','fontsize',16)

subplot(2,2,3)
theta = deg2rad(0:1:360);

g = reshape(Ax, [size(x1),3]);
y = reshape(xAx, size(x1));
% contour(x1,x2,y)
hold
idx = 1:4:size(x1,1);
quiver(x1(idx,idx),x2(idx,idx),-g(idx,idx,1),...
    -g(idx,idx,2))
set(gca,'PlotBoxAspectRatio',[1 1 1])
title('$Ax$','interpreter','latex','fontsize',16)
plot(cos(theta),sin(theta), 'linewidth', 2)
subplot(2,2,4)

g = reshape(xAxAx, [size(x1),3]);
y = reshape(xAxxAx, size(x1));
% contour(x1,x2,y)
hold
idx = 1:2:size(x1,1);
quiver(x1(idx,idx),x2(idx,idx),-g(idx,idx,1),...
    -g(idx,idx,2))
set(gca,'PlotBoxAspectRatio',[1 1 1])
title('$4x^TAx^TAx$','interpreter','latex','fontsize',16)
plot(cos(theta),sin(theta), 'linewidth', 2)