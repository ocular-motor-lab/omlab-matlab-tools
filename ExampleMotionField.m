%%


w = [0 0 0]';
v = [1 0.5 0]';
eyeel = 0;
h =  1;
N = 250;

stim = "Ground plane"; % "Ground plane" or "Sphere at 1m"

[motionField, visualDirections, motionFieldLinear, motionFieldRotational, Jv, Jw] = Geometry3D.CalculateMotionField(N, w, v, h, eyeel, stim);

% plot jacobians (templates)
Jacobians = {Jv,Jw};
x = visualDirections(:,1);
y = visualDirections(:,2);
z = visualDirections(:,3);
dim = {'x','y','z'};
figure
for i=1:2
    for j=1:3

        subplot(2,3,i*3-2+j-1,'nextplot','add')
        axis equal;
        set(gca,'xlim',[-1 1], 'ylim',[-1 1])
        quiver(y.*acosd(x)/90,z.*acosd(x)/90, squeeze(Jacobians{i}(1,j,:)), squeeze(Jacobians{i}(2,j,:)),'linewidth',2)
        if ( i==1)
            xlabel( ['\partial u / \partial ' dim{j} ],'fontsize',14);
            ylabel( ['\partial v / \partial ' dim{j} ],'fontsize',14);
        else
            xlabel( ['\partial u / \omega_' dim{j} ' dt'],'fontsize',14);
            ylabel( ['\partial v / \omega_' dim{j} ' dt'],'fontsize',14);
        end
        set(gca,'xtick',[],'ytick',[])

    end
end

% plot motion fields
figure
axes('nextplot' , 'add')
colors = get(gca,'colororder');
hs = [];
hs.quivertJwflat = quiver(0,0, 0*2 ,0*2,'color',colors(3,:),'linewidth',1,'LineStyle',':');
hs.quivertJvflat = quiver(0,0, 0*2 ,0*2,'color',colors(4,:),'linewidth',1,'LineStyle',':');
hs.quivertAllvflat  = quiver(0,0, 0*2 ,0*2,'color','k','linewidth',1);
legend([hs.quivertJvflat hs.quivertJwflat hs.quivertAllvflat],{'Linear motion' 'Rotational motion' 'Total motion'})


az = y.*acos(x);
el = z.*acos(x);
azdeg = rad2deg(az);
eldeg = rad2deg(el);
motionFieldRotational = rad2deg(motionFieldRotational);
motionFieldLinear = rad2deg(motionFieldLinear);
motionField = rad2deg(motionField);
motionFieldRotational(abs(azdeg)>80 | abs(eldeg)>80) = nan;
motionFieldLinear(abs(azdeg)>80 | abs(eldeg)>80) = nan;
motionField(abs(azdeg)>80 | abs(eldeg)>80) = nan;

set(hs.quivertJwflat, 'xdata', azdeg(:), 'ydata', eldeg(:), 'UData', motionFieldRotational(:,1), 'VData', motionFieldRotational(:,2));
set(hs.quivertJvflat, 'xdata', azdeg(:), 'ydata', eldeg(:), 'UData', motionFieldLinear(:,1), 'VData', motionFieldLinear(:,2));
set(hs.quivertAllvflat, 'xdata', azdeg(:), 'ydata', eldeg(:), 'UData', motionField(:,1), 'VData', motionField(:,2));