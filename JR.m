% Figure to illustrate different retinal coordinate systems.
% https://www.jennyreadresearch.com/code/ReadPhillipsonGlennerster09/

clear all
close all

fs=14; % fontsize to use for axis labels
fsa=18; % fontsize to use for panel labels
fst = 18; % fontsize to use for titles

% Head-centered coordinates of viewed object
objectP = [-6 7 10];

% Coordinates of real eye:
interoculardistance = 6;
eyes.radius = 1.75;
eyes.eyeballcolor = [0.75 0.75 0.5];
eyes.eyeballtransparency = 0.8;
% Specify eye positions
eyes.eyepos.nodal = [0; 0 ; 0];
eyes.eyepos.azimuth = 0;
eyes.eyepos.elevation = 0;
eyes.eyepos.torsion = 0;
% Colors in which to draw coordinate systems
eyes.azcoordcol = 'k' ;
eyes.elcoordcol = 'k' ;

% Specify the coordinate system to be painted on the retina
coordsystem.nlines = 13;
spacing = 15;
coordsystem.sca  = 1.01; % amount by which to exaggerate eyeball radius so coord lines are visible
coordsystem.linewidth = 0.5;
coordsystem.linewidththick = 1.5; % used for plotting meridians
coords = {'CartesianVirtualPlane' 'azlongellong' 'azlongellat' 'azlatellong' 'azlatellat'};
names = {'(x,y)' '(\alpha, \eta)' '(\alpha, \kappa)' '(\beta, \eta)' '(\beta, \kappa)' };
figure('pos',[   5         329        1160         292])
% Draw the eyes
DrawEyes(eyes);
for jcoord=1:length(coords)    
    if jcoord==1
        coordsystem.spacing = eyes.radius * tand(spacing);
        coordsystem.retcol = [0.3 0.3 0.3];
    else
        coordsystem.spacing = spacing;
    end
    subplot(1,length(coords),jcoord)
    % Mark on coordinate system:
    coordsystem.name = coords{jcoord};
    DrawEyes(eyes);
    eyecoordhandles = DrawEyeCoords(eyes,coordsystem);

    camup([0 1 0])
    camtarget([ 0 0 0])
    campos([-7 7 -25])
    camva(7)
    lighting phong
    camlight headlight
    camlight right
    axis equal tight off
    xlim([-1 1]*3.5)
    ylim([-1 1]*3.5)
    zlim([-1 1]*3.5)
    set(gcf,'color','w')
    text(0,1.8*eyes.radius,-eyes.radius,names{jcoord},'fontsize',14)
    text(eyes.radius,1.8*eyes.radius,-eyes.radius,char(64+jcoord),'fontsize',18,'fontw','bold')
    
end
print -dtiff -r800 -zbuffer Fig_DiffRetCoords


function eye=DrawEyes(eyes)
% r = radius of eyeballs

ebr = eyes.radius;
if isfield(eyes,'nsurface')
    [x z y] = sphere(eyes.nsurface);
else
    [x z y] = sphere(50);
end
x = x.*ebr; y=y.*ebr; z = z.*ebr;
for thiseye=1:length(eyes.eyepos)
    [f,v,c]=surf2patch(x + eyes.eyepos(thiseye).nodal(1),y + eyes.eyepos(thiseye).nodal(2),z + eyes.eyepos(thiseye).nodal(3),z);
    eye(thiseye) = patch('faces',f,'vertices',v,'facecol',eyes.eyeballcolor,'edgecol','none','facealpha',eyes.eyeballtransparency,'edgealpha',eyes.eyeballtransparency);
    % If asked, draw in an iris and pupil
    if isfield(eyes,'iriscolor')
        X = get(eye(thiseye),'Xdata')-eyes.eyepos(thiseye).nodal(1);
        Y = get(eye(thiseye),'Ydata')-eyes.eyepos(thiseye).nodal(2);
        Z = get(eye(thiseye),'Zdata')-eyes.eyepos(thiseye).nodal(3);
        [n,m]=size(Z);
        c = repmat(eyes.eyeballcolor,m,1);
        Rot = inv(RotationMatrix(eyes.eyepos(thiseye)));
        R = Rot * [X(1,:);Y(1,:);Z(1,:)];
        Z = R(3,:);
        j = find(Z>eyes.irisradius);
        c(j,:) = repmat(eyes.iriscolor,length(j),1);
        j = find(Z>eyes.pupilradius);
        c(j,:) = repmat(eyes.pupilcolor,length(j),1);
        set(eye(thiseye),'facecolor','flat','edgecol','none','facevertexcdata',c);
    end
end
end



function eyecoordhandles = DrawEyeCoords(eyes,coordsystem)
% Draw an appropriate coordinate system on the spherical eyeballs.

LEFT=1;RIGHT=2;
h=[];
for thiseye=1:length(eyes.eyepos)
    Rot{thiseye}=RotationMatrix(eyes.eyepos(thiseye));
end
% Fields of coordsystem:
%    .nlines = number of gridlines
%    .spacing = spacing in between each line
n = coordsystem.nlines;
gridlinepos = ([1:n]-(n+1)/2)*coordsystem.spacing;
gridline = [ min(gridlinepos) : range(gridlinepos)/100 : max(gridlinepos) ];
switch coordsystem.name
    case 'Cartesian'
        for k=1:2 % First draw lines of varying x, constant y; Next draw lines of varying y, constant x;
            if length(gridline>0)
                for j=1:n
                    % First get positions on cyclopean eye centered at origin, in primary position:
                    nrm = sqrt(gridline.^2 + gridlinepos(j)^2 + eyes.radius^2 );
                    if k==1
                        % First draw lines of varying x, constant y:
                        X = eyes.radius * gridline./ nrm;
                        Y = eyes.radius * gridlinepos(j)./ nrm;
                        col = eyes.azcoordcol;
                    else
                        % Next draw lines of varying y, constant x:
                        Y = eyes.radius * gridline./ nrm;
                        X = eyes.radius * gridlinepos(j)./ nrm;
                        col = eyes.elcoordcol;
                    end
                    Z = -eyes.radius^2./ nrm;
                    if sind(gridlinepos(j))==0
                        lw = coordsystem.linewidththick;
                    else
                        lw = coordsystem.linewidth;
                    end
                    % Now convert to actual eye position:
                    for thiseye=1:length(eyes.eyepos)
                        tmp = coordsystem.sca * Rot{thiseye} * [X;Y;Z] + repmat(eyes.eyepos(thiseye).nodal,1,length(X));
                        hold on
                        h = [h plot3(tmp(1,:),tmp(2,:),tmp(3,:),'col',col,'linew',lw)];
                    end
                end
            end
        end

    case 'CartesianVirtualPlane'
        for k=1:2 % First draw lines of varying x, constant y; Next draw lines of varying y, constant x;
            if length(gridline>0)
                for j=1:n
                    % First get positions on cyclopean eye centered at origin, in primary position:
                    nrm = sqrt(gridline.^2 + gridlinepos(j)^2 + eyes.radius^2 );
                    if k==1
                        % First draw lines of varying x, constant y:
                        X = eyes.radius * gridline./ nrm;
                        Y = eyes.radius * gridlinepos(j)./ nrm;
                        col = eyes.azcoordcol;
                    else
                        % Next draw lines of varying y, constant x:
                        Y = eyes.radius * gridline./ nrm;
                        X = eyes.radius * gridlinepos(j)./ nrm;
                        col = eyes.elcoordcol;
                    end
                    if sind(gridlinepos(j))==0
                        lw = coordsystem.linewidththick;
                    else
                        lw = coordsystem.linewidth;
                    end
                    % Plot lines on virtual plane
                    hold on
                    Z = -eyes.radius^2./ nrm;
                    lam = -eyes.radius./Z;
                    plot3(lam.*X,lam.*Y,-eyes.radius*ones(size(gridline)),'col',col,'linew',lw)
                    % Now convert to actual eye position:
                    for thiseye=1:length(eyes.eyepos)
                        tmp = coordsystem.sca * Rot{thiseye} * [X;Y;Z] + repmat(eyes.eyepos(thiseye).nodal,1,length(X));
                        hold on
                        h = [h plot3(tmp(1,:),tmp(2,:),tmp(3,:),'col',coordsystem.retcol,'linew',lw)];
                    end
                end
            end
        end

    case 'azlongellong'
        for k=1:2 % First draw lines of varying x, constant y; Next draw lines of varying y, constant x;
            if length(gridline>0)
                for j=1:n
                    % First get positions on cyclopean eye centered at origin, in primary position:
                    nrm = sqrt(tand(gridline).^2 + tand(gridlinepos(j))^2 + 1 );
                    if k==1
                        % First draw lines of varying x, constant y:
                        X = eyes.radius * tand(gridline)./ nrm;
                        Y = eyes.radius * tand(gridlinepos(j))./ nrm;
                        col = eyes.azcoordcol ;
                    else
                        % Next draw lines of varying y, constant x:
                        Y = eyes.radius * tand(gridline)./ nrm;
                        X = eyes.radius * tand(gridlinepos(j))./ nrm;
                        col = eyes.elcoordcol ;
                    end
                    Z = -eyes.radius./ nrm * sign(cosd(gridlinepos(j)));
                    if sind(gridlinepos(j))==0
                        lw = coordsystem.linewidththick;
                    else
                        lw = coordsystem.linewidth;
                    end
                    % Now convert to actual eye position:
                    for thiseye=1:length(eyes.eyepos)
                        tmp = coordsystem.sca * Rot{thiseye} * [X;Y;Z] + repmat(eyes.eyepos(thiseye).nodal,1,length(X));
                        hold on
                        h = [h plot3(tmp(1,:),tmp(2,:),tmp(3,:),'col',col,'linew',lw)];
                    end
                end
            end
        end
    case 'azlongellat'
        for k=1:2 % First draw lines of varying x, constant y; Next draw lines of varying y, constant x;
            if length(gridline>0)
                for j=1:n
                    % First get positions on cyclopean eye centered at origin, in primary position:
                    if k==1
                        % First draw lines of varying az, constant el:
                        az = gridline;
                        el = gridlinepos(j)*ones(1,length(gridline));
                        col = eyes.azcoordcol;
                    else
                        % Next draw lines of varying el, constant az:
                        az = gridlinepos(j)*ones(1,length(gridline));
                        el = gridline;
                        col = eyes.elcoordcol;
                    end
                    X = eyes.radius * cosd(el) .* sind(az);
                    Y = eyes.radius * sind(el) ;
                    Z = -eyes.radius * cosd(el) .* cosd(az);
                    if abs(sind(gridlinepos(j)))<1e-10
                        lw = coordsystem.linewidththick;
                    else
                        lw = coordsystem.linewidth;
                    end
                    % Now convert to actual eye position:
                    for thiseye=1:length(eyes.eyepos)
                        tmp = coordsystem.sca * Rot{thiseye} * [X;Y;Z] + repmat(eyes.eyepos(thiseye).nodal,1,length(X));
                        hold on
                        h = [h plot3(tmp(1,:),tmp(2,:),tmp(3,:),'col',col,'linew',lw)];
                    end
                end
            end
        end
    case 'azlatellong'
        for k=1:2 % First draw lines of varying x, constant y; Next draw lines of varying y, constant x;
            if length(gridline>0)
                for j=1:n
                    % First get positions on cyclopean eye centered at origin, in primary position:
                    if k==1
                        % First draw lines of varying az, constant el:
                        az = gridline;
                        el = gridlinepos(j)*ones(1,length(gridline));
                        col = eyes.azcoordcol;
                    else
                        % Next draw lines of varying el, constant az:
                        az = gridlinepos(j)*ones(1,length(gridline));
                        el = gridline;
                        col = eyes.elcoordcol;
                    end
                    X = eyes.radius .* sind(az);
                    Y = eyes.radius * cosd(az).*sind(el) ;
                    Z = -eyes.radius * cosd(az) .* cosd(el);                    
                    if abs(sind(gridlinepos(j)))<1e-10
                        lw = coordsystem.linewidththick;
                    else
                        lw = coordsystem.linewidth;
                    end
                    % Now convert to actual eye position:
                    for thiseye=1:length(eyes.eyepos)
                        tmp = coordsystem.sca * Rot{thiseye} * [X;Y;Z] + repmat(eyes.eyepos(thiseye).nodal,1,length(X));
                        hold on
                        h = [h plot3(tmp(1,:),tmp(2,:),tmp(3,:),'col',col,'linew',lw)];
                    end
                end
            end
        end
    case 'azlatellat'
        for k=1:2 % First draw lines of varying x, constant y; Next draw lines of varying y, constant x;
            if length(gridline>0)
                for j=1:n
                    % First get positions on cyclopean eye centered at origin, in primary position:
                    if k==1
                        % First draw lines of varying az, constant el:
                        az = gridline;
                        el = gridlinepos(j)*ones(1,length(gridline));
                        col = eyes.azcoordcol;
                    else
                        % Next draw lines of varying el, constant az:
                        az = gridlinepos(j)*ones(1,length(gridline));
                        el = gridline;
                        col = eyes.elcoordcol;
                    end
                    X = eyes.radius .* sind(az);
                    Y = eyes.radius .* sind(el) ;
                    Z = -eyes.radius * sqrt(cosd(az).^2 - sind(el).^2);        
                    if abs(sind(gridlinepos(j)))<1e-10
                        lw = coordsystem.linewidththick;
                    else
                        lw = coordsystem.linewidth;
                    end
                    % Remove bits that go out of range
                    nrm = X.^2+Y.^2+Z.^2;
                    jjj = find(nrm<=eyes.radius^2 + 1e-6 & abs(X)<=eyes.radius & abs(Y)<=eyes.radius & abs(Z)<eyes.radius & imag(Z)==0);
                    X=X(jjj);
                    Y=Y(jjj);
                    Z=Z(jjj);
                    % Now convert to actual eye position:
                    for thiseye=1:length(eyes.eyepos)
                        tmp = coordsystem.sca * Rot{thiseye} * [X;Y;Z] + repmat(eyes.eyepos(thiseye).nodal,1,length(X));
                        hold on
                        h = [h plot3(tmp(1,:),tmp(2,:),tmp(3,:),'col',col,'linew',lw)];
                    end
                end
            end
        end


end

% Mark fovea
if isfield(coordsystem,'markfovea')
    for thiseye=1:length(eyes.eyepos)
        tmp = coordsystem.sca * Rot{thiseye} * [0;0;-eyes.radius] + eyes.eyepos(thiseye).nodal;
        h = [h plot3(tmp(1,:),tmp(2,:),tmp(3,:),'col',eyes.azcoordcol,'marker','o','markerfacecol',eyes.azcoordcol,'markersize',12)];
    end
end
if isfield(coordsystem,'markpupil')
    for thiseye=1:length(eyes.eyepos)
        tmp = coordsystem.sca * Rot{thiseye} * [0;0;eyes.radius] + eyes.eyepos(thiseye).nodal;
        h = [h plot3(tmp(1,:),tmp(2,:),tmp(3,:),'col',eyes.azcoordcol,'marker','o','markerfacecol',eyes.azcoordcol,'markersize',12)];
    end
end

eyecoordhandles = h;

end


function Rot=RotationMatrix(eyepos)
% eyepos can either be a single eyeposition structure, or a cell array of
% such structures

if ~iscell(eyepos)
	tmp = eyepos;
	clear eyepos;
	eyepos{1} = tmp;
end

for j=1:length(eyepos)
	cosH = cosd(eyepos{j}.azimuth);
	sinH = sind(eyepos{j}.azimuth);
	cosV = cosd(eyepos{j}.elevation);
	sinV = sind(eyepos{j}.elevation);
	cosT = cosd(eyepos{j}.torsion);
	sinT = sind(eyepos{j}.torsion);

	% Find rotation matrices
	RotH = [cosH 0 sinH ; 0 1 0 ; -sinH 0 cosH];
	RotV = [ 1 0 0 ; 0 cosV -sinV ; 0  sinV cosV ];
	RotT = [cosT -sinT 0 ; sinT cosT 0; 0  0 1];
	Rot{j} = RotV*RotH*RotT;
end

if exist('tmp','var') % then the original eyepos supplied was not a cell array; return something that is not a cell array
	Rot = Rot{1};
end
end