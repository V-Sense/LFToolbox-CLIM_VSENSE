% LFLytroDecodeImage - decode a Lytro light field from a raw lenslet image, called by LFUtilDecodeLytroFolder
%
% Usage:
%     [LF, LFMetadata, WhiteImageMetadata, LensletGridModel, DecodeOptions] = ...
%         LFLytroDecodeImage( InputFname, DecodeOptions )
%     [LF, LFMetadata, WhiteImageMetadata, LensletGridModel, DecodeOptions] = ...
%         LFLytroDecodeImage( InputFname )
%
% This function decodes a raw lenslet image into a 4D light field. Its purpose is to tailor the core lenslet decoding
% function, LFDecodeLensletImageSimple, for use with Lytro data. It is envisioned that other camera formats will be
% supported by similar functions in future.
% 
% Supported file formats include Lytro LFP files and extracted .raw files accompanied by metadata, as extracted by using
% LFP Reader v2.0.0, for example. See LFToolbox.pdf for more information.
%
% The white image appropriate to a light field is selected based on a white image database, and so
% LFUtilProcessWhiteImages must be run before this function.
% 
% LFUtilDecodeLytroFolder is useful for decoding multiple images.
%
% The optional DecodeOptions argument includes several fields defining filename patterns. These are
% combined with the input LFFnameBase to build complete filenames. Filename patterns include the
% placeholder '%s' to signify the position of the base filename. For example, the defualt filename
% pattern for a raw input file, LensletImageFnamePattern, is '%s__frame.raw'. So a call of the
% form LFLytroDecodeImage('IMG_0001') will look for the raw input file 'IMG_0001__frame.raw'.
% 
% Inputs:
% 
%     InputFname :  Filename of the input light field -- the extension is used to detect LFP or raw input.
% 
%   [optinal] DecodeOptions : all fields are optional, defaults are for LFP Reader v2.0.0 naming
%         .WhiteProcDataFnameExtension : Grid model from LFUtilProcessWhiteImages, default 'grid.json'
%          .WhiteRawDataFnameExtension : White image file extension, default '.RAW'
%              .WhiteImageDatabasePath : White image database, default 'Cameras/WhiteImageDatabase.mat'
% 
%            For compatibility with extracted .raw and .json files:
%                .MetadataFnamePattern : JSON file containing light field metadata, default '_metadata.json'
%              .SerialdataFnamePattern : JSON file containing serial numbers, default '_private_metadata.json'
%
% Outputs:
% 
%                     LF : 5D array containing a 4-channel (RGB + Weight) light field, indexed in
%                          the order [j,i,l,k, channel]
%             LFMetadata : Contents of the metadata and serial metadata files
%     WhiteImageMetadata : Conents of the white image metadata file
%      LensletGridModel  : Lenslet grid model used to decode the light field, as constructed from
%                          the white image by LFUtilProcessWhiteImages / LFBuildLensletGridModel
%          DecodeOptions : The options as applied, including any default values omitted in the input
% 
% 
% Example:
%     
%     LF = LFLytroDecodeImage('Images/F01/IMG_0001__frame.raw');
%     or
%     LF = LFLytroDecodeImage('Images/Illum/LorikeetHiding.lfp');
%     
%     Run from the top level of the light field samples will decode the respective raw or lfp light fields.
%     LFUtilProcessWhiteImages must be run before decoding will work.
% 
% See also: LFUtilDecodeLytroFolder, LFUtilProcessWhiteImages, LFDecodeLensletImageSimple, LFSelectFromDatabase

% Part of LF Toolbox v0.4 released 12-Feb-2015
% Copyright (c) 2013-2015 Donald G. Dansereau

% modified by Rodrigo Daudt
%        - modifications of white balancing and gamma parameters in the case of ILLUM cameras to be consistent with F01 cameras.
% modified by Mikael Le Pendu : 25 Aug. 2016
%        - added loading of black images (used for hot pixel correction).
%        - New decoded option HotPixelCorrect.

function [LF, LFMetadata, WhiteImageMetadata, LensletGridModel, DecodeOptions] = ...
    LFLytroDecodeImage( InputFname, DecodeOptions )

%---Defaults---
DecodeOptions = LFDefaultField( 'DecodeOptions', 'WhiteProcDataFnameExtension', '.grid.json' );
DecodeOptions = LFDefaultField( 'DecodeOptions', 'WhiteRawDataFnameExtension', '.RAW' );
DecodeOptions = LFDefaultField( 'DecodeOptions', 'WhiteImageDatabasePath', fullfile('Cameras','WhiteImageDatabase.mat'));
% Compatibility: for loading extracted raw / json files
DecodeOptions = LFDefaultField( 'DecodeOptions', 'MetadataFnamePattern', '_metadata.json' );
DecodeOptions = LFDefaultField( 'DecodeOptions', 'SerialdataFnamePattern', '_private_metadata.json' );
DecodeOptions = LFDefaultField( 'DecodeOptions', 'HotPixelCorrect', true );

%---
LF = [];
LFMetadata = [];
WhiteImageMetadata = [];
LensletGridModel = [];

%---Read the LFP or raw file + external metadata---
FileExtension = InputFname(end-2:end);
switch( lower(FileExtension) )
    case 'raw' %---Load raw light field and metadata---
        FNameBase = InputFname(1:end-4);
        MetadataFname = [FNameBase, DecodeOptions.MetadataFnamePattern];
        SerialdataFname = [FNameBase, DecodeOptions.SerialdataFnamePattern];
        fprintf('Loading lenslet image and metadata:\n\t%s\n', InputFname);
        LFMetadata = LFReadMetadata(MetadataFname);
        LFMetadata.SerialData = LFReadMetadata(SerialdataFname);
        
        switch( LFMetadata.camera.model )
            case 'F01'
                BitPacking = '12bit';
                DecodeOptions.DemosaicOrder = 'bggr';
            case 'B01'
                BitPacking = '10bit';
                DecodeOptions.DemosaicOrder = 'grbg';
        end
        LensletImage = LFReadRaw(InputFname, BitPacking);
        
    otherwise %---Load Lytro LFP format---
        fprintf('Loading LFP %s\n', InputFname );
        LFP = LFReadLFP( InputFname );
        if( ~isfield(LFP, 'RawImg') )
            fprintf('No light field image found, skipping...\n');
            return
        end
        LFMetadata = LFP.Metadata;
        LFMetadata.SerialData = LFP.Serials;
        LensletImage = LFP.RawImg;
        DecodeOptions.DemosaicOrder = LFP.DemosaicOrder;
end

%---Select appropriate white image---
DesiredCam = struct('CamSerial', LFMetadata.SerialData.camera.serialNumber, ...
    'ZoomStep', LFMetadata.devices.lens.zoomStep, ...
    'FocusStep', LFMetadata.devices.lens.focusStep );
DecodeOptions.WhiteImageInfo = LFSelectFromDatabase( DesiredCam, DecodeOptions.WhiteImageDatabasePath );
PathToDatabase = fileparts( DecodeOptions.WhiteImageDatabasePath );
if( isempty(DecodeOptions.WhiteImageInfo) || ~strcmp(DecodeOptions.WhiteImageInfo.CamSerial, DesiredCam.CamSerial) )
    fprintf('No appropriate white image found, skipping...\n');
    return
end

%---Display image info---
%---Check serial number---
fprintf('\nWhite image / LF Picture:\n');
fprintf('%s, %s\n', DecodeOptions.WhiteImageInfo.Fname, InputFname);
fprintf('Serial:\t%s\t%s\n', DecodeOptions.WhiteImageInfo.CamSerial, DesiredCam.CamSerial);
fprintf('Zoom:\t%d\t\t%d\n', DecodeOptions.WhiteImageInfo.ZoomStep, DesiredCam.ZoomStep);
fprintf('Focus:\t%d\t\t%d\n\n', DecodeOptions.WhiteImageInfo.FocusStep, DesiredCam.FocusStep);

%---Load white image, white image metadata, and lenslet grid parameters---
DecodeOptions.WhiteImageInfo.Fname = strrep(DecodeOptions.WhiteImageInfo.Fname, '\', '/'); % FIX FOR LINUX
WhiteMetadataFname = fullfile(PathToDatabase, DecodeOptions.WhiteImageInfo.Fname);
WhiteProcFname = LFFindLytroPartnerFile(WhiteMetadataFname, DecodeOptions.WhiteProcDataFnameExtension);
WhiteRawFname = LFFindLytroPartnerFile(WhiteMetadataFname, DecodeOptions.WhiteRawDataFnameExtension);

fprintf('Loading white image and metadata...\n');
LensletGridModel = LFStruct2Var(LFReadMetadata(WhiteProcFname), 'LensletGridModel');
WhiteImageMetadataWhole = LFReadMetadata( WhiteMetadataFname );
WhiteImageMetadata = WhiteImageMetadataWhole.master.picture.frameArray.frame.metadata;
WhiteImageMetadata.SerialData = WhiteImageMetadataWhole.master.picture.frameArray.frame.privateMetadata;

switch( WhiteImageMetadata.camera.model )
    case 'F01'
        assert( WhiteImageMetadata.image.rawDetails.pixelPacking.bitsPerPixel == 12 );
        assert( strcmp(WhiteImageMetadata.image.rawDetails.pixelPacking.endianness, 'big') );
        DecodeOptions.LevelLimits = [LFMetadata.image.rawDetails.pixelFormat.black.gr, LFMetadata.image.rawDetails.pixelFormat.white.gr];
        DecodeOptions.ColourMatrix = reshape(LFMetadata.image.color.ccmRgbToSrgbArray, 3,3);
        DecodeOptions.ColourBalance = [...
            LFMetadata.image.color.whiteBalanceGain.r, ...
            LFMetadata.image.color.whiteBalanceGain.gb, ...
            LFMetadata.image.color.whiteBalanceGain.b ];
        
        if(DecodeOptions.CLIMLegacyCols)
            DecodeOptions.SensorNormalizeRBGains = [1 1];
            DecodeOptions.Gamma = LFMetadata.image.color.gamma^0.5;
        else
            %DecodeOptions.SensorNormalizeRBGains =   [LFMetadata.devices.sensor.normalizedResponses.gr/LFMetadata.devices.sensor.normalizedResponses.r, ...
            %                                          LFMetadata.devices.sensor.normalizedResponses.gr/LFMetadata.devices.sensor.normalizedResponses.b];
            DecodeOptions.SensorNormalizeRBGains = [1 1];
            DecodeOptions.Gamma = 'sRGB';
        end
        
        DecodeOptions.ExposureBias = LFMetadata.image.modulationExposureBias + 1;
        
        BitPacking = '12bit';
        
    case 'B01'
        assert( WhiteImageMetadata.image.rawDetails.pixelPacking.bitsPerPixel == 10 );
        assert( strcmp(WhiteImageMetadata.image.rawDetails.pixelPacking.endianness, 'little') );
        DecodeOptions.LevelLimits = [LFMetadata.image.pixelFormat.black.gr, LFMetadata.image.pixelFormat.white.gr];
        DecodeOptions.ColourBalance = ([...
            LFMetadata.image.color.whiteBalanceGain.r, ...
            LFMetadata.image.color.whiteBalanceGain.gb, ...
            LFMetadata.image.color.whiteBalanceGain.b ]);
        DecodeOptions.SatLevels = [WhiteImageMetadata.devices.sensor.normalizedResponses.r, WhiteImageMetadata.devices.sensor.normalizedResponses.gr, WhiteImageMetadata.devices.sensor.normalizedResponses.b];
        %Normalized gains (The White Image should have the same multipliers)
        if(DecodeOptions.CLIMLegacyCols)
            DecodeOptions.SensorNormalizeRBGains = [1 1];
            DecodeOptions.Gamma = 1/2.4^0.5; %R.DAUDT consistency with gamma for F01 cameras.
            DecodeOptions.ColourMatrix = (reshape(LFMetadata.image.color.ccm, 3,3) + 1.0*eye(3))/2.0; %Added by R.DAUDT : 1.0*eye(3))/2.0
%            DecodeOptions.ColourMatrix = reshape(LFMetadata.image.color.ccm, 3,3);
            DecodeOptions.ColourBalance = sqrt(DecodeOptions.ColourBalance);
        else
            DecodeOptions.SensorNormalizeRBGains =   [WhiteImageMetadata.devices.sensor.normalizedResponses.gr/WhiteImageMetadata.devices.sensor.normalizedResponses.r, ...
                                                     WhiteImageMetadata.devices.sensor.normalizedResponses.gr/WhiteImageMetadata.devices.sensor.normalizedResponses.b];
            DecodeOptions.Gamma = 'sRGB';
            %Conversion Matrice from Non-normalized sensor data to sRGB.
            DecodeOptions.ColourMatrix = reshape(WhiteImageMetadata.image.color.ccmRgbToSrgbArray, 3,3);
        end
        
        DecodeOptions.ExposureBias = LFMetadata.image.modulationExposureBias;
        
        BitPacking = '10bit';
        
    otherwise
        fprintf('Unrecognized camera model, skipping...\');
        return
end
WhiteImage = LFReadRaw( WhiteRawFname, BitPacking );

if(DecodeOptions.HotPixelCorrect)
    %---Select and read black images---
    % Assumes that the black image numbers are 0 and 1 and that they are
    % located in the same folder as the white images.
    [ImagesDir, BlackImagesNames] = LFSearchBlackImages(WhiteRawFname,[0 1]);
    BlackImageSum = 0;
    for i=1:length(BlackImagesNames)
        BlackImageSum = BlackImageSum + double(LFReadRaw( [ImagesDir '\' BlackImagesNames{i}], BitPacking ));
    end

    %---Compute Hot pixels---
    m = mean(BlackImageSum(:));
    HotPixels = false(size(WhiteImage));
    HotPixels(BlackImageSum > 2*m) = 1;
    clear BlackImageSum
else
    HotPixels = [];
end


%---Decode---
fprintf('Decoding lenslet image :');
[LF, LFWeight, DecodeOptions] = LFDecodeLensletImageSimple( LensletImage, WhiteImage, HotPixels, LensletGridModel, DecodeOptions );
LF(:,:,:,:,end+1:end+DecodeOptions.NWeightChans) = LFWeight;
