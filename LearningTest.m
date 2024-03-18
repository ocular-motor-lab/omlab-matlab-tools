%%
dt = 0.01;
t = (0:dt:15)';

close all


n = 100;
m = 2;
A=eye(m+1);
A(end) = -1;

S0 = randn(n,m);


x = zeros(length(t),m);
x(t>=0 & t<3,:) = randn(size(x(t>=0 & t<3,:)));
x(t>=3 & t<5,1) = 1;
x(t>=5 & t<7,2) = 1;
x(t>=7,:) = randn(size(x(t>=7,:)));
x = x./repmat(vecnorm(x,2,2),1,m);

[t, sout] = ode45(@(ti,xi)AttractorNetworkUpdateS(ti, xi, interp1(t,x,ti)', A), t, S0(:)); 


Sout = nan(n,m,length(t));
for i=1:length(t)
    Sout(:,:,i) = reshape(sout(i,:),n,m);
end

figure

subplot(3,5,1:5);
plot(t, x,'Linewidth',2); ylim([-2 2]), title('input')

subplot(3,5,6:10);
plot(t, sout); title('coefficients')

t2 = find(t>3, 1,'first');
t3= find(t>5, 1,'first');
t4 = find(t>7, 1,'first');

subplot(3,5,11);
plot(Sout(:,1,1), Sout(:,2,1),'o'), title('t1 start')
subplot(3,5,12);
plot(Sout(:,1,t2), Sout(:,2,t2),'o'), title('t2 end of random')
subplot(3,5,13);
plot(Sout(:,1,t3), Sout(:,2,t3),'o'), title('t3 end of first bias')
subplot(3,5,14);
plot(Sout(:,1,t4), Sout(:,2,t4),'o'), title('t4 end of second bias')
subplot(3,5,15);
plot(Sout(:,1,end), Sout(:,2,end),'o'), title('t5 end')

h = get(gcf,'children');
set(h(1:5),'xlim',[-1 1]*1.5, 'ylim',[-1 1]*1.5)