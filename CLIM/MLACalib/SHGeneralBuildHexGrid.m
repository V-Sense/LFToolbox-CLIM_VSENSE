% July 7 2016 put this in a separate file so I can run
% SHGeneralBuildHexGrid.m
% A general function that can build a hex grid
% operating on either the toolbox LensletGridModel or on the
% augmented LensletGridModel with the centers added in.
% 
% SHGeneralBuildHexGrid.m uses .CenterX and .CenterY fields in
% LensletGridModel if they exist to return the calibrated grid.  This
% LensletGridModel corresponds to that stored in the .gridcent.json
% files, and uses the (centerx, centery) coordinate as the rotation
% point.
%
% If a LensletGridModel is missing this information, the ...
% function builds a hex grid that is rotated around the origin.  
% Note that this grid is what is returned using LFBuildHexGrid if 
% called with the LensletGridModel in the original metadata files
% .grid.json, OR with the calibrated LensletGridModel (without
% center information) stored in the metadata files .grid2.json.
%
% If .grid2.json metadata is used, the returned coordinates are NOT
% exactly aligned with those of the perfectly calibrated grid (due to
% the difference in the rotation point, origin v. center respectively).

function [GridCoords] = SHGeneralBuildHexGrid( LensletGridModel )

% if there is no center in the LensletGridModel
if (myIsField(LensletGridModel, 'CenterX')) == 0

  %sprintf('No center specified; will rotate around origin')
  RotCent = eye(3);
  RotCent(1:2,3) = [LensletGridModel.HOffset, LensletGridModel.VOffset];

  ToOffset = eye(3);
  ToOffset(1:2,3) = [LensletGridModel.HOffset, LensletGridModel.VOffset];

  R = ToOffset * RotCent * LFRotz(LensletGridModel.Rot) * RotCent^-1;


  [vv,uu] = ndgrid((0:LensletGridModel.VMax-1).*LensletGridModel.VSpacing, (0:LensletGridModel.UMax-1).*LensletGridModel.HSpacing);

  uu(LensletGridModel.FirstPosShiftRow:2:end,:) = uu(LensletGridModel.FirstPosShiftRow:2:end,:) + 0.5.*LensletGridModel.HSpacing;

  GridCoords = [uu(:), vv(:), ones(numel(vv),1)];
  GridCoords = (R*GridCoords')';

  GridCoords = reshape(GridCoords(:,1:2), [LensletGridModel.VMax,LensletGridModel.UMax,2]);

else

  % According to SHLFBuildCalibGrid.m - the HSpacing and VSpacing
  % that are in the jsoncent file are on the square grid.  

  sprintf('Center specified; will rotate around center')
  Wdim = 3280;
  Htim = 3280;

  cx = LensletGridModel.CenterX;
  cy = LensletGridModel.CenterY;
  xlength = LensletGridModel.HSpacing;
  ydistYrot = LensletGridModel.VSpacing;
  rotangle = LensletGridModel.Rot;

  xnum = floor(cx/xlength);
  ynum = floor(cy/ydistYrot);

  xdirvect = cx-xnum*xlength:xlength:Wdim;
  ydirvect = cy-ynum*ydistYrot: ydistYrot:Htim;
  [vind,uind] = ndgrid(ydirvect, xdirvect); 

  gridcenterx = find(abs(xdirvect - cx)<10^-9);
  gridcentery = find(abs(ydirvect - cy)<10^-9);
  if 2*floor(gridcentery/2) == gridcentery
    shiftstart = 1;
  else
    shiftstart = 2;
  end

  uind(shiftstart:2:end,:) = uind(shiftstart:2:end,:) + xlength/2;

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

  GridCoords = LytroGridCoords;

end





