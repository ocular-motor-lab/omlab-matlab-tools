v = VideoReader('IMG_1076.MOV'); % 6dof
%v = VideoReader('C:\Users\jom\Downloads\IMG_1076.MOV'); % 6dof
% v = VideoReader('C:\Users\jom\Downloads\IMG_1072.MOV'); % Torsion
% v = VideoReader('C:\Users\jom\Downloads\IMG_1071.MOV');
% v = VideoReader('C:\Users\jom\Downloads\IMG_1070.MOV');
close all
[Xo,Yo] = meshgrid(1:v.Width,1:v.Height);
X = Xo(:);Y=Yo(:);

xc = mean(Xo(:));
yc = mean(Yo(:));

Xo = Xo - xc;
Yo = Yo - yc ;


I1 = double(im2gray(v.readFrame()));
I2 = double(im2gray(v.readFrame()));


% spherical coordinates of pixels
sensorWidth = 4; %mm
% focalLengthEquiv = 13; %mm
% cropFactor = 8.6;
hRes = v.Width;
AOVh = 108;
d = (sensorWidth/2)/tand(AOVh/2);


if (1)
    %%
    STD = 300;
    N = 500;

    p = exp(-(Xo.^2 + Yo.^2)/(STD^2));
    r = rand(size(Xo));
    l = log(r./p);
    [s,i] = sort(l(:));

    isel = i(1:N);

    X = Xo(isel);
    Y = Yo(isel);
    Isel = I1(isel);
    figure

    % F1 = scatteredInterpolant(X(:), Y(:), I1(:));
    % Isel = F1(Xsel, Ysel);
    F = scatteredInterpolant(X, Y, Isel,'nearest');
    Iinterp = F(Xo(:),Yo(:));
    Iinterp = reshape(Iinterp, size(I1));

    figure
    subplot(2,1,1)
    imagesc(I1)
    % subplot(3,1,2)
    % DT = delaunayTriangulation(Xsel,Ysel);
    % voronoi(Xsel,Ysel)
    % hold
    % triplot(DT)
    subplot(2,1,2)
    imagesc(Iinterp)
    colormap gray

    % hold
    % scatter(X(:)+xc,Y(:)+yc,'r.')
    % scatter(Xdif+xc,Ydif+yc,'g.')
    
end

    DT = delaunayTriangulation(X,Y);
    DTedges = DT.edges;
    Xdif = (X(DTedges(:,1))+X(DTedges(:,2)))/2;
    Ydif = (Y(DTedges(:,1))+Y(DTedges(:,2)))/2;

% right handed spherical coordinates with x pointing forward
% S is points on the sphere in cartesian coordinates
% Sa is in polar (Angles)
S = [ repmat(d,numel(X),1) -X(:)/hRes*sensorWidth Y(:)/hRes*sensorWidth];
S = S./repmat(vecnorm(S,2,2),1,3);
Sa = [acosd(S(:,2)) acosd(S(:,3))];

% bidirectional connections
L1 = [DTedges(:,1);DTedges(:,2)];
L2 = [DTedges(:,2);DTedges(:,1)];
% angular distance
L12ah = Sa(L2,1) - Sa(L1,1);
L12av = Sa(L2,2) - Sa(L1,2);
L12ax = S(L2,1) - S(L1,1);
L12ay = S(L2,2) - S(L1,2);
L12az = S(L2,3) - S(L1,3);



RV = cross(S,repmat([1 0 0],height(S),1));
Q = quaternion(RV, 'rotvec');
R = Q.rotmat("point");


% LS = [[DTedges(:,1);DTedges(:,2)], [DTedges(:,2);DTedges(:,1)], );
return


% I Need to be able to get two jacobians
%
% One is the jacobian of pixels over angles
%
% Another is the jacobian of pixels over linear displacements
%
% I might as well just generalize a bit more and go directly to 
% random locations of pixels. The math will be exactly the same and cannot
% approximate anything useful really. 
%
% The only thing that can be approximated is the 2D foveal motion

% gradients of the angles to convert the gradients to degs
% this should mean something like what is the torsion, horizontal, or
% vertifcal rotation in degrees that would move you from pixel to pixel
% [Gtx,Gty] = imgradientxy(reshape(S(:,1),size(I1)));

 [Gxx,Gxy] = imgradientxy(acosd(reshape(S(:,2),size(I1))),"central");
 [Gyx,Gyy] = imgradientxy(acosd(reshape(S(:,3),size(I1))),"central");
 [Gtx] = atan2d(Gxx,Gxy);
  Gty = atan2d(Gxy,Gyy);

 [GtxL,GtyL] = imgradientxy((reshape(S(:,1),size(I1))),"central");
 [GxxL,GxyL] = imgradientxy((reshape(S(:,2),size(I1))),"central");
 [GyxL,GyyL] = imgradientxy((reshape(S(:,3),size(I1))),"central");

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




angvel = nan(v.NumFrames,3);
transvel = nan(v.NumFrames,3);
rmsI = nan(v.NumFrames,3);
vout = VideoWriter('Test6dof.mp4', 'MPEG-4');

open(vout);


    Rx = squeeze(R(3,:,:))';
    y = squeeze(R(2,:,:))';

    Tx = squeeze(R(2,:,:))';
    Ty = squeeze(R(3,:,:))';

Rx = [Gtx(:) ,Gxx(:) , Gyx(:) ];
Ry = [Gty(:) ,Gxy(:) , Gyy(:) ];
Tx = [GtxL(:) ,GxxL(:) , GyxL(:) ]/100;
Ty = [GtyL(:) ,GxyL(:) , GyyL(:) ]/100;

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
            quiver3(S(i,1),S(i,2),S(i,3),Gtx(i)/2,Gxx(i)/2,Gyx(i)/2, 'linewidth',2,'color','r')
            % translational x axis % rotational y axis
            quiver3(S(i,1),S(i,2),S(i,3),Gty(i)/2,Gxy(i)/2,Gyy(i)/2, 'linewidth',2,'color','b')
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
h(2).Parent.CLim = [-100 100]*0.1;
h(3).Parent.CLim = [-100 100]*0.1;
h(4).Parent.CLim = [-100 100]*0.1;
h(5).Parent.CLim = [-100 100]*0.1;
h(6).Parent.CLim = [-100 100]*0.1;



i=0;
t2 = v.CurrentTime;
while hasFrame(v)
    i = i+1;

    I1 = double(im2gray(readFrame(v)));
    I1 = imgaussfilt(I1,3);
    t1 = v.CurrentTime;
    [Gx,Gy] = imgradientxy((I1+I2)/2,"central");
    Gx = imgaussfilt(Gx,3);
    Gy = imgaussfilt(Gy,3);
    dI = (I1-I2);
    dt = (t1-t2);
    dIdt = dI/dt;
    t2 = t1;



    % convert the gradients to espherical coordinates
     GxR = Gx;%./Gxx;
     GyR = Gy;%./Gyy;
     GxL = Gx;%./GxxL/250;
     GyL = Gy;%./GyyL/250;


    GR = GxR(:).*Rx + GyR(:).*y;
    GL = GxL(:).*Tx + GyL(:).*Ty;

    [angvel(i,:), ~,Res] = regress(dIdt(:),GR);
    I2 = I1;

    ResI = reshape(Res, size(I1));

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
