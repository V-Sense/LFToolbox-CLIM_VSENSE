function [mlaData] = load_mla_calib(WIDir)
% Load the calibration data from Lytro for the camera with direcory WIDir.

mlaFileNameF01 = 'T1CALIB__MLACALIBRATION.TXT';
mlaFileNameIllum = 'mla_calibration.json';

if(exist([WIDir '/' mlaFileNameF01],'file'))
    mlaFile = mlaFileNameF01;
    Cameratype = 'F01';
elseif(exist([WIDir '/' mlaFileNameIllum],'file'))
    mlaFile = mlaFileNameIllum;
    Cameratype = 'ILLUM';
else
    mlaData = [];
    return;
end

mlaData = LFReadMetadata( [WIDir '/' mlaFile] );
mlaData.zoomPosition = mlaData.zoomPosition';
mlaData.numFocusPositions = mlaData.numFocusPositions';

numZoom = length(mlaData.zoomPosition);
numFocusMax = max(mlaData.numFocusPositions);
fields = fieldnames(mlaData);
for i=1:numel(fields)
    if (~strcmp(fields{i},'zoomPosition') && ~strcmp(fields{i},'numFocusPositions') )
        if (iscell(mlaData.(fields{i})))
            mlaData.(fields{i}) = unknownToNan(cellTo2DMat(mlaData.(fields{i}),[numZoom,numFocusMax]), mlaData.numFocusPositions);
        elseif (isnumeric(mlaData.(fields{i})))
            mlaData.(fields{i}) = unknownToNan(mlaData.(fields{i}), mlaData.numFocusPositions);
        end
    end
end
