%% Figures for appendix on random walks


%% Figure 1 - simulating simple white noise
% https://www.gaussianwaves.com/2013/11/simulation-and-analysis-of-white-noise-in-matlab/

clear all; clc; close all;
Fs = 100;
Ts = 1/Fs;
L = 20000;   % Sample length for the random signal
t = (0:L-1)*Ts;
mu = 0;       % Mean of the random signal  
sigma = 2;    % Standard deviation of the random signal
X = sigma*randn(L,1);

figure();
subplot(2,2,1)
plot(t,X);
title(['White noise : \mu_x=',num2str(mu),' \sigma^2=',num2str(sigma^2), ' Fs=', num2str(Fs)])
xlabel('Time (s)')
ylabel('Sample Values of x(n)')
grid on;


subplot(2,2,2)
n=30; %number of Histrogram bins
[f,x]=hist(X,n);
bar(x,f/trapz(x,f),'EdgeColor','none'); hold on;
%Theoretical PDF of Gaussian Random Variable
g=(1/(sqrt(2*pi)*sigma))*exp(-((x-mu).^2)/(2*sigma^2));
plot(x,g,'linewidth',2); grid on;
title('Theoretical PDF and Simulated Histogram of White Gaussian Noise');
legend('Histogram','Theoretical PDF');
xlabel('x');
ylabel('PDF f_x(x)');


subplot(2,2,3,'nextplot','add')
Rxx=1/L*conv(flipud(X),X);
lags=(-L+1):1:(L-1);

%Alternative method
%[Rxx,lags] =xcorr(X,'biased'); 
%The argument 'biased' is used for proper scaling by 1/L
%Normalize auto-correlation with sample length for proper scaling

plot(lags,Rxx); 
h = plot(0,sigma.^2,'ro');
legend(h,{'Theoretical R(0)'});
title('Auto-correlation Function of white noise');
xlabel('Lags')
ylabel('Correlation')
grid on;


subplot(2,2,4,'nextplot','add')
%Verifying the constant PSD of White Gaussian Noise Process
%with arbitrary mean and standard deviation sigma

N = 1024; 
% split the signal into windows to average power across each of them
z = reshape(X(1:(N*floor(length(X)/N))),N,floor(length(X)/N))';

%By default, FFT is done across each column - Normal command fft(z)
%Finding the FFT across each row
Z = 1/sqrt(N*Fs)*fft(z,[],2); %Scaling by sqrt(N) to get psd
Pzavg = mean(Z.*conj(Z),1);%Computing the mean power from fft
Pzavg=fftshift(Pzavg); %Shift zero-frequency component to center of spectrum
freq=(-N/2:N/2-1)/N*Fs;

% pwelch(X,[],[],N,[],'centered') would give a similar answer

h1 = plot(freq,(Pzavg));
h2 = plot(freq,sigma.^2/Fs*ones(size(freq)),'linewidth',2);
legend([h1 h2],{'Average periodogram', 'Theoretical PSD' })
% plot(normFreq,10*log10(Pzavg),'r');
axis([-0.5*Fs 0.5*Fs 0 10/Fs]); grid on;
ylabel('Power Spectral Density (linear)');
xlabel('Frequency (1/s)');
title('Power spectral density of white noise');

%% Figure 2 - Wiener process

r = round(rand(100,1))*2-1;
clear x
x(1000*(1:100)) = r;
x = cumsum(x);

figure('color','w')
subplot(1,2,1);
plot((1:length(x))/1000,x);
xlabel('Time step')
ylabel('x[n]')
title('Random walk')

x = round(rand(30*1000,1))*2-1;
x = cumsum(x/1000);

subplot(1,2,2);
plot((1:length(x))/1000,x);
xlabel('Time (s)')
ylabel('x(t)')
title('Wiener process')

%% Figure 3 - Spectrum (not perfect)
close all
f = 0:0.001:100;
w = 2*pi*f;
t1 = 0.02;
t2 = 20;
D = 1;

x = 2*D ./ (1/t2^2 + w.^2 +t1^2*w.^4);
v = w.^2*2*D ./ (1/t2^2 + w.^2 +t1^2*w.^4);

figure('color','w')
subplot(1,2,1)
loglog(f,x)
hold

yl = get(gca,'ylim');
line([1 1]*1/(t1*t2), yl)
line([1 1]*1/(t2), yl)
xlabel('Frequency (Hz)')
ylabel('S_x (\omega)')
text(1/(t1*t2), yl(2), '1/\tau_2\tau_1','VerticalAlignment','top')
text((1/t2), yl(2), '1/\tau_2','VerticalAlignment','top')

subplot(1,2,2)
loglog(f,v)
hold

yl = get(gca,'ylim');
line([1 1]*1/(t1*t2), yl)
line([1 1]*1/(t2), yl)
xlabel('Frequency (Hz)')
ylabel('S_v (\omega)')

text(1/(t1*t2), yl(2), '1/\tau_2\tau_1','VerticalAlignment','top')
text((1/t2), yl(2), '1/\tau_2','VerticalAlignment','top')

%% 0 - Simulations of simple random walks (Wiener process)
clear all, close all
Reps = 5000;

Ts = 1/100;    % sampling rate
D = 40;         % diffusion constant in arcmin^2/s
VarV = 2*D/Ts;
SigmaV = sqrt(VarV);
SigmaVdeg = SigmaV/60;

t = (0:Ts:2)';

v = SigmaV * randn(height(t),Reps);
x = cumsum(v*Ts);
% Note this is is equivalent to x = filter(Ts,[1 -1], v);



VarX = std(x,[],2).^2;
vmeasured = diff(x)./diff(t);


figure('color','w')
subplot(2,2,1);
plot(t,x(:,1:min(20,size(x,2))));
xlabel('Time (s)')
ylabel('Position (arcmin)')
title(sprintf('Brownian motion (Fs = %0.1f Hz)',1/Ts))


subplot(2,2,2);
plot(t,VarX,'o');
hold
% plot(t,VarXmeasured,'o');
plot(t,2*D*t,'r','linewidth',2)
set(gca,'xlim',[0 max(t)/2])
xlabel('Time (s)')
ylabel('MSD (arcmin^2)')
title(sprintf('Mean square displacements (D = %0.1f minarc^2)',D))

subplot(2,2,3);
binwidth = 4/60;
vBins = -10:binwidth:10;
vBinsC = (vBins(1:end-1)+vBins(2:end))/2;
bar(vBinsC, histcounts(v(:)/60,vBins)/length(v(:)),'EdgeColor','none');
hold
% bar(vBinsC, histcounts(vmeasured(:),vBins)/length(vmeasured(:)));
plot(vBinsC,sqrt(1/(2*pi*SigmaVdeg^2))*exp(-vBinsC.^2/(2*SigmaVdeg^2))*binwidth,'r','linewidth',2)
title('Velocity distribution')
xlabel('Velocity (deg/s)')
ylabel('Probability')

N = 128;
X = vmeasured(:);
z = reshape(X(1:(N*floor(length(X)/N))),N,floor(length(X)/N))';

%By default, FFT is done across each column - Normal command fft(z)
%Finding the FFT across each row
Z = 1/sqrt(N/Ts)*fft(z,[],2); %Scaling by sqrt(N) to get psd
Pzavg = mean(Z.*conj(Z),1);%Computing the mean power from fft
Pzavg=fftshift(Pzavg); %Shift zero-frequency component to center of spectrum
freq=(-N/2:N/2-1)/N*1/Ts;

subplot(2,2,4)
bar(freq, Pzavg,'EdgeColor','none');
hold
plot(freq,2*D*ones(size(freq)),'r','linewidth',2)

 set(gca,'ylim',[0 2*60])
title('Velocity power spectrum')
xlabel('Frequency (Hz)')
ylabel('Velocity PSD (arcmin^2/s/Hz)')
legend({'Simulation' 'Theoretical'},'box','off')


 % 1 - Simulations of simple random walks with relaxation time cosntant

%% varying taus
% clear all, close all

figure
tiledlayout(4,2,'padding', "tight", "TileSpacing", "tight");

Reps = 5000;
dur =2; % duration of each rep
Fs = 1000;
D = 40;         % diffusion constant arcmin^2/s

tau1 = 0.0005; % Zuber
tau2 = 20; % not sure
[x, t] = GenerateRandomWalk1D(Fs, dur, D, tau1, tau2, Reps);
PlotRW(t, x, D, tau1, tau2,Fs)


tau1 = 0.05; % Zuber
tau2 = 20; % not sure
[x, t] = GenerateRandomWalk1D(Fs, dur, D, tau1, tau2, Reps);
PlotRW(t, x, D, tau1, tau2,Fs)


tau1 = 0.5; % Zuber
tau2 = inf; % not sure
[x, t] = GenerateRandomWalk1D(Fs, dur, D, tau1, tau2, Reps);
PlotRW(t, x, D, tau1, tau2,Fs)


tau1 = 0.05; % Zuber
tau2 = 0.5; % not sure
[x, t] = GenerateRandomWalk1D(Fs, dur, D, tau1, tau2, Reps);
x = x + randn(size(x));
PlotRW(t, x, D, tau1, tau2,Fs)

legend({'Simulation','2Dt','Short timescale aprox','Long timescale steady state'})





%% adding noise and bias
% clear all, close all

figure
tiledlayout(3,2,'padding', "tight", "TileSpacing", "tight");

tau1 = 0.05; % Zuber
tau2 = 1; % not sure
[x, t] = GenerateRandomWalk1D(Fs, dur, D, tau1, tau2, Reps);
PlotRW(t, x, D, tau1, tau2,Fs)

x1 = x + randn(size(x))*2;
PlotRW(t, x1, D, tau1, tau2,Fs)

bias = 0.01;
x2 = x + cumsum(bias*ones(size(x)));
PlotRW(t, x2, D, tau1, tau2,Fs)

legend({'Simulation','2Dt','Short timescale aprox','Long timescale steady state'})


function PlotRW(t, x, D, tau1, tau2,Fs)


VarX = mean(x.^2,2);
vmeasured = diff(x)./diff(t);
Ts = 1/Fs;
SigmaV = sqrt(2*D/Ts);


SteadyStateVarx = D*tau2;


nexttile
plot(t,x(:,1:min(20,size(x,2))));
xlabel('Time (s)')
ylabel('Position (arcmin)')
title(sprintf('Brownian motion (Fs = %0.1f Hz, {\\tau_1} = %0.0fms, {\\tau_2} = %0.1f s)',1/Ts,tau1*1000, tau2))


%     2*D * (t + tau1*( exp(-t/tau1) - 1 ) ),...
nexttile
plot(t,VarX,'o');
hold
% plot(t,VarXmeasured,'o');
plot(t,2*D*t,'r--','linewidth',1)
plot(t,...
    2*D*tau1*(  t/tau1 - 3/2  +2*exp(-t/tau1)  -1/2*exp(-2*t/tau1)  ), ...
    'r','linewidth',2)
if ( SteadyStateVarx< D*3)
    plot(t,SteadyStateVarx*ones(size(t)),'g','linewidth',2)
end

set(gca,'xlim',[0 max(t)/2])
xlabel('Time (s)')
ylabel('MSD (ardmin^2)')
title(sprintf('Mean square displacements (D=%0.1f minarc^2)',D))

return
nexttile
vBins = -150:150;
vBinsC = (vBins(1:end-1)+vBins(2:end))/2;
bar(vBinsC, histcounts(vmeasured(:),vBins)/length(vmeasured(:)),'EdgeColor','none')
hold
% bar(vBinsC, histcounts(vmeasured(:),vBins)/length(vmeasured(:)));
plot(vBinsC, ...
    sqrt(1/(2*pi*SigmaV^2))*exp(-vBinsC.^2/(2*SigmaV^2)) ...
    , 'r--','linewidth',1)
SigmaVt = sqrt(D/tau1); % only works if tau >> Ts;
plot(vBinsC,sqrt(1/(2*pi*SigmaVt^2))*exp(-vBinsC.^2/(2*SigmaVt^2)),'r','linewidth',2)
title('Velocity distribution ~N(0, 2DF_s)')
xlabel('Velocity (deg/s)')
ylabel('Probability')


N = 1024;
X = vmeasured(:);
z = reshape(X(1:(N*floor(length(X)/N))),N,floor(length(X)/N))';

%By default, FFT is done across each column - Normal command fft(z)
%Finding the FFT across each row
Z = 1/sqrt(N/Ts)*fft(z,[],2); %Scaling by sqrt(N) to get psd
Pzavg = mean(Z.*conj(Z),1);%Computing the mean power from fft
Pzavg=fftshift(Pzavg); %Shift zero-frequency component to center of spectrum
freq=(-N/2:N/2-1)/N*1/Ts;

nexttile
bar(freq, Pzavg,'EdgeColor','none')
hold
w = 2*pi*freq;
plot(freq, 2*D./(tau1^2*w.^2+1),'r','linewidth',2)
plot(freq, 2*D*ones(size(freq)),'r--','linewidth',1)
set(gca,'ylim',[0 2])
title('Velocity power spectrum {= 2D / (\tau^2w^2+1)}')
xlabel('Frequency (Hz)')
ylabel('Velocity PSD (deg^2/s/Hz)')

legend({'Simulation' 'Theoretical' 'Theoretical for no {\tau}'},'box','off')

end

function [x, t] = GenerateRandomWalk1D(Fs, t, D, tau1, tau2, Reps)

if ( nargin == 0 )
    t = 2;          % duration
    Fs = 100;       % sampling rate
    tau1 = 0.046;   % plant time constant - Zuber:46ms
    tau2 = 20;      % integrator time constant - 20s
    D = 40;         % diffusion constant arcmin^2/s
end

if (~exist('Reps','var'))
    Reps = 1; % repeatitions
end


Ts = 1/Fs;

VarV = 2*D/Ts;
SigmaV = sqrt(VarV);

t = (0:Ts:t)';
v = SigmaV*randn(numel(t),Reps);

b = Ts^2;
a = [(tau1 + Ts + Ts^2/tau2) -(2*tau1+Ts) +tau1];
x = filter(b,a,v);

% equivalent to filter(b,a,v)
% x = zeros(numel(t),Reps);
% for i=1:length(t)
%     if ( i>3 )
%         x(i,:) = 1/a(1)*(b(1)*v(i,:) -a(2)*x(i-1,:) -a(3)*x(i-2,:) );
%     elseif ( i>2 )
%         x(i,:) = 1/a(1)*(b(1)*v(i,:) -a(2)*x(i-1,:) );
%     else
%         x(i,:) = 1/a(1)*(b(1)*v(i,:) );
%     end
% end

end

