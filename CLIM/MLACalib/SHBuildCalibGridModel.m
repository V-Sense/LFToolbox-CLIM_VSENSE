% SHBuildCalibGridModel.m
% July 5 2016
% build the structure CalibGridModel for a given white image
%
% Retrieves mlaCalib data corresponding to a given zoom and focus.
% Only returns data required by SHLFBuildCalibGrid.m to compute the lenslet grid model.
% The data is added to the output structure CalibGridModel:
% - CalibGridModel.xDiskOriginPixels
% - CalibGridModel.yDiskOriginPixels
% - CalibGridModel.xDiskStepPixelsX
% - CalibGridModel.xDiskStepPixelsY
% - CalibGridModel.yDiskStepPixelsX
% - CalibGridModel.yDiskStepPixelsY

function [CalibGridModel] = SHBuildCalibGridModel(ZoomStep, FocusStep, mla_data)
if(isempty(mla_data))
    CalibGridModel=[];
else
    zoomId = find(mla_data.zoomPosition==ZoomStep, 1);
    focusId = find(mla_data.focusPosition(zoomId,:)==FocusStep, 1);
    if(isempty(zoomId) || isempty(focusId))
        CalibGridModel=[];
    else
        CalibGridModel.xDiskOriginPixels = mla_data.xDiskOriginPixels(zoomId,focusId);
        CalibGridModel.yDiskOriginPixels = mla_data.yDiskOriginPixels(zoomId,focusId);
        CalibGridModel.xDiskStepPixelsX = mla_data.xDiskStepPixelsX(zoomId,focusId);
        CalibGridModel.xDiskStepPixelsY = mla_data.xDiskStepPixelsY(zoomId,focusId);
        CalibGridModel.yDiskStepPixelsX = mla_data.yDiskStepPixelsX(zoomId,focusId);
        CalibGridModel.yDiskStepPixelsY = mla_data.yDiskStepPixelsY(zoomId,focusId);
    end
end