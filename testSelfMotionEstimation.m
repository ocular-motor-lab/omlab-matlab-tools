v = VideoReader('IMG_1076.MOV'); % 6dof
%v = VideoReader('C:\Users\jom\Downloads\IMG_1076.MOV'); % 6dof
% v = VideoReader('C:\Users\jom\Downloads\IMG_1072.MOV'); % Torsion
% v = VideoReader('C:\Users\jom\Downloads\IMG_1071.MOV');
% v = VideoReader('C:\Users\jom\Downloads\IMG_1070.MOV');
close all
[X,Y] = meshgrid(1:v.Width,1:v.Height);
X = X - mean(X(:));
Y = Y - mean(Y(:));


% spherical coordinates of pixels
sensorWidth = 4; %mm
% focalLengthEquiv = 13; %mm
% cropFactor = 8.6;
hRes = v.Width;
AOVh = 108;
d = (sensorWidth/2)/tand(AOVh/2);
% right handed spherical coordinates with x pointing forward
S = [ repmat(d,numel(X),1) -X(:)/hRes*sensorWidth Y(:)/hRes*sensorWidth];
S = S./repmat(vecnorm(S,2,2),1,3);


% gradients of the angles to convert the gradients to degs
% this should mean something like what is the torsion, horizontal, or
% vertifcal rotation in degrees that would move you from pixel to pixel
% [Gtx,Gty] = imgradientxy(reshape(S(:,1),size(I1)));
 [Gxmag] = imgradient(reshape(S(:,2),size(I1)));
 [Gymag] = imgradient(reshape(S(:,3),size(I1)));
 [Gxx,Gxy] = imgradientxy(acosd(reshape(S(:,2),size(I1))),"central");
 [Gyx,Gyy] = imgradientxy(acosd(reshape(S(:,3),size(I1))),"central");

% %%
% figure
% subplot(3,3,1)
% imagesc(reshape(S(:,1),size(I1'))');
% subplot(3,3,2)
% imagesc(Gtx);
% subplot(3,3,3)
% imagesc(Gty);
% 
% subplot(3,3,4)
% imagesc(reshape(S(:,2),size(I1'))');
% subplot(3,3,5)
% imagesc(Gxx);
% subplot(3,3,6)
% imagesc(Gxy);
% 
% subplot(3,3,7)
% imagesc(reshape(S(:,2),size(I1'))');
% subplot(3,3,8)
% imagesc(Gyx);
% subplot(3,3,9)
% imagesc(Gyy);




RV = cross(S,repmat([1 0 0],height(S),1));
Q = quaternion(RV, 'rotvec');
R = Q.rotmat("point");

angvel = nan(v.NumFrames,3);
transvel = nan(v.NumFrames,3);
rmsI = nan(v.NumFrames,3);
vout = VideoWriter('Test6dof.mp4');

open(vout);


if ( 0)
    %%
    % figure

    figure
    hold
    r = randi(height(S),100);
    for ir=1:length(r)
         i=r(ir);
            % vector pointing out
            quiver3(S(i,1),S(i,2),S(i,3),R(1,1,i)/10,R(1,2,i)/10,R(1,3,i)/10, 'linewidth',2,'color','r')
            % translational x axis % rotational y axis
            quiver3(S(i,1),S(i,2),S(i,3),R(2,1,i)/10,R(2,2,i)/10,R(2,3,i)/10, 'linewidth',2,'color','b')
            % rotational x axis % translational y axis
            quiver3(S(i,1),S(i,2),S(i,3),R(3,1,i)/10,R(3,2,i)/10,R(3,3,i)/10, 'linewidth',2,'color','g')
    end
        view(45,30)
    set(gca,'xlim',[-1 1])
    set(gca,'ylim',[-1 1])
    set(gca,'zlim',[-1 1])
    set(gca,'PlotBoxAspectRatio',[1 1 1])

     %%
    % figure
close all
    figure
    hold
    r = randi(height(S),200);
    for ir=1:length(r)
         i=r(ir);
            % vector pointing out
            quiver3(S(i,1),S(i,2),S(i,3),Gtx(i)/10,Gxx(i)/10,Gyx(i)/10, 'linewidth',2,'color','r')
            % translational x axis % rotational y axis
            quiver3(S(i,1),S(i,2),S(i,3),Gty(i)/10,Gxy(i)/10,Gyy(i)/10, 'linewidth',2,'color','b')
            % % rotational x axis % translational y axis
            % quiver3(S(i,1),S(i,2),S(i,3),R(3,1,i)/10,R(3,2,i)/10,R(3,3,i)/10, 'linewidth',2,'color','g')
    end
        view(45,30)
    set(gca,'xlim',[-1 1])
    set(gca,'ylim',[-1 1])
    set(gca,'zlim',[-1 1])
    set(gca,'PlotBoxAspectRatio',[1 1 1])
end


%%

I1 = double(im2gray(v.readFrame()));
I2 = double(im2gray(v.readFrame()));

clear h;

f=figure('color','w');

pos = f.Position;
f.Position = [pos(1)-(1040-pos(3)) pos(2)-(780-pos(4))  1040 780];

subplot(4,3,1)
h(1) = imshow(I1);
title('Current frame')
subplot(4,3,4)
dIdt = I2-I1;
h(2) = imagesc(dIdt);
colormap(redblue)
title('Change over time')

[Gx,Gy] = imgradientxy(I1);
subplot(4,3,2)
h(3) = imagesc(Gx);
title('Horizontal gradient')
subplot(4,3,3)
h(4) = imagesc(Gy);
title('Vertical gradient')


subplot(4,3,5)
dIdt = I2-I1;
h(5) = imagesc(dIdt);
colormap(redblue)
title('Res. change over time after rot')

subplot(4,3,6)
dIdt = I2-I1;
h(6) = imagesc(dIdt);
colormap(redblue)
title('Res. change over time after lin')


subplot(4,3,[7 8 9])
hl = plot(angvel);
ylabel('Ang Velocity');
legend({'H' 'V' 'T'})
ylim([-1 1]*50);
xlim([0 v.NumFrames])

subplot(4,3,[10 11 12])
hl2 = plot(transvel);
ylabel('Lin Velocity');
legend({'Fw' 'lat' 'vert'})
ylim([-1 1]*50);
xlim([0 v.NumFrames])

h(1).Parent.CLim = [0 255];
h(2).Parent.CLim = [-100 100];
h(3).Parent.CLim = [-100 100];
h(4).Parent.CLim = [-100 100];
h(5).Parent.CLim = [-100 100];
h(6).Parent.CLim = [-100 100];



    Rx = squeeze(R(3,:,:))';
    Ry = squeeze(R(2,:,:))';

    Tx = squeeze(R(2,:,:))';
    Ty = squeeze(R(3,:,:))';

i=0;
t2 = v.CurrentTime;
while hasFrame(v)
    i = i+1;

    I1 = double(im2gray(readFrame(v)));
    I1 = imgaussfilt(I1,5);
    t1 = v.CurrentTime;
    [Gx,Gy] = imgradientxy(I1,"central");
    Gx = imgaussfilt(Gx,5);
    Gy = imgaussfilt(Gy,5);
    dI = (I1-I2);
    dt = (t1-t2);
    dIdt = dI/dt;
    t2 = t1;



    % convert the gradients to espherical coordinates
     GxR = Gx./Gxx;
     GyR = Gy./Gyy;


    GR = GxR(:).*Rx + GyR(:).*Ry;
    [angvel(i,:), ~,Res] = regress(dIdt(:),GR);
    I2 = I1;

    ResI = reshape(Res, size(I1));

    GL = Gx(:).*Tx + Gy(:).*Ty;
    [transvel(i,:), ~,Res2] = regress(Res,GL);
    ResI2 = reshape(Res2, size(I1));

    rmsI(i,1) = rms(dIdt(:));
    rmsI(i,2) = rms(Res);
    rmsI(i,3) = rms(Res2);

    set(h(1),'CData',I1);
    set(h(2),'CData',dI);
    set(h(3),'CData',Gx);
    set(h(4),'CData',Gy);
    set(h(5),'CData',ResI*dt);
    set(h(6),'CData',ResI2*dt);
    set(hl(1),'ydata',angvel(:,3));
    set(hl(2),'ydata',angvel(:,2));
    set(hl(3),'ydata',angvel(:,1));
    set(hl2(1),'ydata',transvel(:,1));
    set(hl2(2),'ydata',transvel(:,2));
    set(hl2(3),'ydata',transvel(:,3));
    drawnow

    if ( mod(i,5)==0)
       currFrame = getframe(gcf);
       writeVideo(vout,currFrame);
    end

end
close(vout);
