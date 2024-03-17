%%
dt = 0.001;
t = (0:dt:10)';

close all


n = 100;
m = 2;
A=eye(m+1);
A(end) = -1;

S0 = randn(n,m);


x = zeros(length(t),m);
x(t>=0 & t<5,:) = randn(size(x(t>=0 & t<5,:)));
x(t>=5 & t<7.5,1) = 1;
x(t>=7.5 ,2) = 1;
x = x./repmat(vecnorm(x,2,2),1,m);

[t, sout] = ode45(@(ti,xi)AttractorNetworkUpdateS(ti, xi, interp1(t,x,ti)', A), t, S0(:)); 


Sout = nan(n,m,length(t));
for i=1:length(t)
    Sout(:,:,i) = reshape(sout(i,:),n,m);
end

figure

subplot(3,5,1:5);
plot(x,'Linewidth',2); ylim([-1 2]), title('input')

subplot(3,5,6:10);
plot(sout); title('coefficients')

subplot(3,5,11);
plot(Sout(:,1,1), Sout(:,2,1),'o'), title('Start')
subplot(3,5,12);
plot(Sout(:,1,round(end/4)), Sout(:,2,round(end/4)),'o'), title('by 1/4')
subplot(3,5,13);
plot(Sout(:,1,round(end*2/4)), Sout(:,2,round(end*2/4)),'o'), title('by 2/4')
subplot(3,5,14);
plot(Sout(:,1,round(end*3/4)), Sout(:,2,round(end*3/4)),'o'), title('by 3/4')
subplot(3,5,15);
plot(Sout(:,1,end), Sout(:,2,end),'o'), title('End')

h = get(gcf,'children')
set(h(1:5),'xlim',[-1 1]*1.5, 'ylim',[-1 1]*1.5)