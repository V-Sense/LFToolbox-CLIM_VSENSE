% July 5 2016 
% SHLFBuildCalibGrid.m
% Build a calibrated grid using the calibration information from
% Lytro for a zoom/focus combination, and create an "augmented"
% CalibratedLensletGridModel that has 2 additional fields, 
% CenterX and CenterY.
% 
% This augmented CalibratedLensletGridModel can be passed to
% SHGeneralBuildHexGrid which will return exactly the calibrated
% model that SHLFBuildCalibGrid returns.  (If a LensletGridModel
% without center information is passed to SHGeneralBuildHexGrid,
% it builds a grid exactly as LFBuildHexGrid does.)
%
% For SHLFBuildCalibGrid, pass in a single structure with the
% following calibration information:
%
% CalibGridModel.xDiskOriginPixels
% CalibGridModel.yDiskOriginPixels
% CalibGridModel.xDiskStepPixelsX
% CalibGridModel.xDiskStepPixelsY
% CalibGridModel.yDiskStepPixelsX
% CalibGridModel.yDiskStepPixelsY


function [CalibLensletGridModel, LytroGridCoords] = SHLFBuildCalibGrid( CalibGridModel, Wdim, Hdim )

% center of the grid in pixels
cx = CalibGridModel.xDiskOriginPixels;
cy = CalibGridModel.yDiskOriginPixels;

% get the two vectors that define the hex grid - these vectors
% point to the closest horizontal and vertical (at approx 60 degree
% angle) lenslets from the center.

xdistX = CalibGridModel.xDiskStepPixelsX;
xdistY = CalibGridModel.xDiskStepPixelsY;
ydistX = CalibGridModel.yDiskStepPixelsX;
ydistY = CalibGridModel.yDiskStepPixelsY;

% Compute the spacing of the lenslets, and the rotation that will 
% make the hex grid parallel with x-axis.  Note that my angle is
% positive, while the angle that is computed by the LFToolbox from
% the white images themselves is negative.  My angle should always
% be passed to the BuildHexGrid as positive.

% spacing of lenslets
try
xlength = sqrt(xdistX^2 + xdistY^2);
catch err
    x=1;
end
ylength = sqrt(ydistX^2 + ydistY^2);

angle_from_horiz = atan(xdistY/xdistX);
rotangle = angle_from_horiz;

% get new y coordinates on the rotated grid

% The new angle of the y vector is the old angle, minus the rotation.
newyangle = atan(ydistY/ydistX) - rotangle;

% Now compute new y-coords on rotated grid so x-vector is parallel with x-axis
ydistYrot = ylength*sin(newyangle);
ydistXrot = ylength*cos(newyangle);


% get boundaries for ndgrid 
xnum = floor(cx/xlength);
ynum = floor(cy/ydistYrot);

% Create grid coordinates on a squared grid, meaning that
% the x-vector pointing to the horizontal adjacent lenslet
% has a 0 y-component.

% A big difference between my grid build and that in the toolbox is
% that the toolbox indices ndgrid from 0:Umax-1 or Vmax-1, scaled
% by the lenslet spacing.  I use the actual grid.

xdirvect = cx-xnum*xlength:xlength:Wdim;
ydirvect = cy-ynum*ydistYrot: ydistYrot:Hdim;

% Get center index.  These should correspond to the
% calibrated center.  Also used to figure out hex grid shifts of
% every other row.
gridcenterx = find(abs(xdirvect - cx)<10^-9);
gridcentery = find(abs(ydirvect - cy)<10^-9);

[vind,uind] = ndgrid(ydirvect, xdirvect); 

% For debugging - yes, these are correct in my development.
% Checked - now I have the grid with the identified centers.
% uind(gridcentery, gridcenterx);
% vind(gridcentery, gridcenterx);


% Now shift uuind by (ydistXrot or xlength?) to do the hexagonal tiling.  Execute
% this shift at the center, which is in the middle.  So I have to keep
% the center in the same place.  If gridcentery is even, then shift
% the odds.  Otherwise shift the evens.


if 2*floor(gridcentery/2) == gridcentery
  shiftstart = 1;
else
  shiftstart = 2;
end

% Note we shift by the projection of the rotated yspacing onto the
% x axis, because we are now on a fully rectangular grid.
%uind(shiftstart:2:end,:) = uind(shiftstart:2:end,:) + ydistXrot;

% I ended up deciding to shift by xlength/2.  
uind(shiftstart:2:end,:) = uind(shiftstart:2:end,:) + xlength/2;

% Now add in an extra column if necessary.  This portion of code was
% developed in ComputeLytroGrids.m.

% Now - if we have shifted such that we can fit an extra column in,
% we need to shift all the rows that can fit the extra column.
pp = size(uind);
newpt = uind(shiftstart,1) - xlength;  % the first point that is
                                      % shifted
uindn = uind;
for k=shiftstart:2:pp(1)
  uindn(k,:) = [newpt uind(k, 1:end-1)];
end

clear uind; uind = uindn;
% this now has changed gridcenterx and gridcentery, so don't use
% them anymore.
clear gridcenterx gridcentery;


% ROTATE AROUND THE CENTER, NOT THE ORIGIN
% So subtract the center, do the rotation, and then add the center
% back in.

RotCent = eye(3);
RotCent(1:2,3) = [cx, cy];  % inverse of this matrix will have the
                            % negative of the center point

RotMat = eye(3);
RotMat(1:2,1:2) = [cos(rotangle) -sin(rotangle); sin(rotangle) ...
	  cos(rotangle)];

bigR = RotCent*RotMat*(RotCent^-1);


LytroGridCoords = [uind(:), vind(:), ones(numel(uind), 1)];
LytroGridCoords = (bigR*LytroGridCoords')';

LytroGridCoords = reshape(LytroGridCoords(:,1:2), [size(uind), 2]);

ss = size(LytroGridCoords);

CalibLensletGridModel.HSpacing = xlength;  % on a square grid
CalibLensletGridModel.VSpacing = ydistYrot; % on a square grid
CalibLensletGridModel.HOffset = LytroGridCoords(2,1,1);
CalibLensletGridModel.VOffset = LytroGridCoords(1,1,2);
CalibLensletGridModel.Rot = rotangle;
CalibLensletGridModel.UMax = ss(2);
CalibLensletGridModel.VMax = ss(1);
CalibLensletGridModel.FirstPosShiftRow = shiftstart;
CalibLensletGridModel.Orientation = 'horz';  % this is horz for
                                            % everything, so hard
                                            % code it

% new information
CalibLensletGridModel.CenterX = cx;
CalibLensletGridModel.CenterY = cy;

