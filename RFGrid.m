
% Tangent grid of the receptive field


 defaultparams.gridSizeHdeg = 2;
 defaultparams.gridSizeVdeg = 1;
 defaultparams.gridSeparationHdeg = 0.1;
 defaultparams.gridSeparationVdeg = 0.1;

 defaultparams.RfAz = 10;
 defaultparams.RfEl = 10;

 defaultparams.eyeAz = -15;
 defaultparams.eyeEl = -15;

 defaultparams.DisplaySizeHdeg = 40;
 defaultparams.DisplaySizeVdeg = 30;



app = InteractiveUI('RF GRID',@(app) (UpdatePlot(app.Values, app.Data)), 0.1);
app.AddSlider('eyeAz',          defaultparams.eyeAz,  [-90 90])
app.AddSlider('eyeEl',          defaultparams.eyeEl,  [-90 90])
app.AddSlider('RfAz',           defaultparams.RfAz,  [-90 90])
app.AddSlider('RfEl',           defaultparams.RfEl,  [-90 90])

app.AddSlider('gridSizeHdeg',    defaultparams.gridSizeHdeg,  [0 90] )
app.AddSlider('gridSizeVdeg',   defaultparams.gridSizeVdeg,  [0 90] )
app.AddSlider('gridSeparationHdeg',   defaultparams.gridSeparationHdeg,  [0 10] )
app.AddSlider('gridSeparationVdeg',   defaultparams.gridSeparationVdeg,  [0 10] )

app.AddSlider('DisplaySizeHdeg',    defaultparams.DisplaySizeHdeg,  [0 90] )
app.AddSlider('DisplaySizeVdeg',    defaultparams.DisplaySizeVdeg,  [0 90] )



app.Data.hf = figure;
app.Data.hax = axes('nextplot','add');
app.Data.hgrid1 = plot(0, 0,'+');
app.Data.hgrid2 = plot(0, 0,'r+');
app.Data.hGaze = plot(0, 0,'o','markersize',10, 'LineWidth',3);

set(gca,'xlim', defaultparams.DisplaySizeHdeg*[-1 1], 'ylim',defaultparams.DisplaySizeVdeg*[-1 1])
grid
axis equal


app.Open();

% UpdatePlot(defaultparams, data)

function UpdatePlot(params, data)
% convert to radians to define the grid as tangent to a unit sphere
Xrf = params.gridSizeHdeg*deg2rad(rem(1,params.gridSeparationHdeg) + -1:params.gridSeparationHdeg:1);
Yrf = params.gridSizeVdeg*deg2rad(rem(1,params.gridSeparationHdeg) + -1:params.gridSeparationHdeg:1);
[Xrfgrid, Yrfgrid] = meshgrid(Xrf, Yrf);

% get the vectors from the origin to the tangent grid and normalized them
% so they are in the unit sphere
R = sqrt(Xrfgrid.^2 + Yrfgrid.^2 + 1);
Xrfgrid = Xrfgrid./R;
Yrfgrid = Yrfgrid./R;
ZrfGrid = 1./R;

% Now we change to world coordinates so z is up
Rfgrid3dvec = [ZrfGrid(:)  Yrfgrid(:) Xrfgrid(:)];
z = zeros(size(Rfgrid3dvec(:,1)));

%%
% coordinates of receptive field position and eye position
% azimuth and elevation define the polar coordinates of the orientation
% by defining the eccentricity as the sqrt(Elevation.^2 + Azimuth.^2)
% and the angle as the atan2(Elevation, -Azimuth)

eyePositionRotMat = Geometry3D.List2Mat([atan2(params.eyeAz, params.eyeEl) deg2rad(sqrt(params.eyeEl.^2 + params.eyeAz.^2)) 0]);
rfPositionRotMat = Geometry3D.List2Mat([atan2(params.RfAz, params.RfEl) deg2rad(sqrt(params.RfEl.^2 + params.RfAz.^2)) 0]);

% Rotate grid to by position and then by eye position
rotatedRFVectors = (eyePositionRotMat*rfPositionRotMat*Rfgrid3dvec')';

RFInScreen =  rotatedRFVectors(:,[3 2])./rotatedRFVectors(:,1);

RFinScreenApprox = [Xrfgrid(:)+deg2rad(params.RfAz+params.eyeAz) Yrfgrid(:)+deg2rad(params.RfEl+params.eyeEl) ];

set(data.hgrid1, 'xdata', rad2deg(RFInScreen(:,1)), 'ydata',rad2deg(RFInScreen(:,2)));
set(data.hgrid2, 'xdata', rad2deg(RFinScreenApprox(:,1)), 'ydata', rad2deg(RFinScreenApprox(:,2)));


set(data.hGaze, 'xdata', params.eyeAz, 'ydata',params.eyeEl);

set(gca,'xlim', params.DisplaySizeHdeg*[-1 1], 'ylim',params.DisplaySizeVdeg*[-1 1])
end