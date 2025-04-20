%% relative angle betwen retinal motion and eye motion distributions of objects that moving vs stable in the world 
% 
doangles = 0;

if ( doangles)
    bins = -180:1:180;
else
    bins = -2:0.05:2;
end

eyeMotionStd = 2;
objectMotionStd = 1;
effCopySTd = 0.5;


objectMotion = randn(10000,2)*objectMotionStd;
eyeMotion = randn(10000,2)*eyeMotionStd;
noise = randn(10000,2)*effCopySTd;
retinalMotion =  -eyeMotion + objectMotion;
eff = noise + eyeMotion;
v = dot(eff', retinalMotion');
v1 = v ./ (vecnorm(eff,2,2)'.*vecnorm(retinalMotion,2,2)');
if ( doangles)
 v1 = acosd(v1);
end


eyeMotion = randn(10000,2)*eyeMotionStd;
noise = randn(10000,2)*effCopySTd;
retinalMotion =  -eyeMotion;
eff = noise + eyeMotion;
v = dot(eff', retinalMotion');
v2 = v ./ (vecnorm(eff,2,2)'.*vecnorm(retinalMotion,2,2)');

if ( doangles)
    v2 = acosd(max(min(real(v2),1),-1));
end



figure
histogram(v1, bins);
hold
histogram(v2, bins);


%% relative amplitude betwen retinal motion and eye motion distributions of objects that moving vs stable in the world 

bins = -2:0.05:20;

eyeMotionStd = 2;
objectMotionStd = 1;
effCopySTd = 0.5;


objectMotion = randn(10000,2)*objectMotionStd;
eyeMotion = randn(10000,2)*eyeMotionStd;
noise = randn(10000,2)*effCopySTd;
retinalMotion =  -eyeMotion + objectMotion;
eff = noise + eyeMotion;
v1 = vecnorm(retinalMotion,2,2) ./ vecnorm(eff,2,2)-1;


eyeMotion = randn(10000,2)*eyeMotionStd;
noise = randn(10000,2)*effCopySTd;
retinalMotion =  -eyeMotion;
eff = noise + eyeMotion;
v2 = vecnorm(retinalMotion,2,2) ./ vecnorm(eff,2,2)-1;




figure
histogram(v1, bins);
hold
histogram(v2, bins);



%% relative amplitude betwen retinal motion and eye motion distributions of objects that moving vs stable in the world 

bins = -2:0.01:2;

eyeMotionStd = 1;
objectMotionStd = 2;
effCopySTd = 0.5;

N = 100000;

objectMotion = randn(N,2)*objectMotionStd;
eyeMotion = randn(N,2)*eyeMotionStd;
noise = randn(N,2)*effCopySTd;
retinalMotion =  -eyeMotion + objectMotion;
eff = noise + eyeMotion;
v1s = vecnorm(retinalMotion,2,2) ./ vecnorm(eff,2,2)-1;
v = dot(eff', retinalMotion');
v1a = v ./ (vecnorm(eff,2,2)'.*vecnorm(retinalMotion,2,2)');

v1 = v1s/5 + v1a';

eyeMotion = randn(N,2)*eyeMotionStd;
noise = randn(N,2)*effCopySTd;
retinalMotion =  -eyeMotion;
eff = noise + eyeMotion;
v2s = vecnorm(retinalMotion,2,2) ./ vecnorm(eff,2,2)-1;
v = dot(eff', retinalMotion');
v2a = v ./ (vecnorm(eff,2,2)'.*vecnorm(retinalMotion,2,2)');


v2 = v2s/5 + v2a';
close all
figure
subplot(3,1,1)
histogram(v1a, bins);
hold
histogram(v2a, bins);

subplot(3,1,2)
histogram(v1s, bins);
hold
histogram(v2s, bins);

subplot(3,1,3)
histogram(v1, bins);
hold
histogram(v2, bins);










%%
objectMotion = randn(10000,2)*2;
eyeMotion = randn(10000,2)*2;

figure
hold
for i=-5:5
eyeMotion = randn(10000,2)*0 + i;
eyeMotionNormalized = eyeMotion;
v = dot(eyeMotion',eyeMotion'-objectMotion');
% figure


histogram(v,100);
end