# LFToolbox-CLIM_VSENSE

This is the Light Field toolbox used for lenslet RAW data demultiplexing in [1].
For more details see [our webpage](https://v-sense.scss.tcd.ie/research/light-fields/a-pipeline-for-lenslet-light-field-quality-enhancement/).
The additional post-processing steps presented in [1] and [6] (recolouring for view consistency and denoising) are not present in this toolbox but are available as separate tools: [Recolouring](https://github.com/V-Sense/LFToolbox_Recolouring_HPR), [Denoising](https://github.com/V-Sense/LFBM5D).

This code is built upon the [Light Field toolbox v0.4-CLIM](https://www.irisa.fr/temics/demos/lightField/CLIM/DataSoftware.html), itself based on the [Light Field toolbox v0.4](https://uk.mathworks.com/matlabcentral/fileexchange/49683-light-field-toolbox-v0-4) by Donald Dansereau [2].

When using this code in your research, please cite the papers [1], [6].

## Modifications from Light Field toolbox v0.4-CLIM:

- Corrected white balance and exposure as detailed in [1],[6] (M. LE PENDU)
	- Always enabled.
	- White image normalization (so that devignetting does not interfere with white balance and exposure).
	- Corrected usage of White balance and exposure bias metadata.
	- Apply standard sRGB gamma correction.

- Highlight Processing as detailed in [6] (M. LE PENDU)
	- Enable/Disable with DecodeOptions.CorrectSaturated (default=true).
	- Process saturated pixels for better highlights (e.g. corrects pink highlights issue and retrieve some details).

- Prevent any clipping in the pipeline (M. LE PENDU)
	- Enable/Disable with DecodeOptions.noClip (default=false).
	- Generates a variable maxLum. The final data in LF must be multiplied by maxLum to recover original dynamic.
	- If disabled, a soft clipping is performed to avoid loosing all the details in the highlights.

- Early White Balance (M. LE PENDU)
	- Enable/Disable with DecodeOptions.EarlyWhiteBalance (default=true).
	- Perform white balance directly on the RAW data, before demosaicing.

- New interpolation mode: 'none' (M. LE PENDU)
	- Enable/Disable with DecodeOptions.ResampMethod='none' (default='fast').
	- Non destructive mode: generates many incomplete views with mask indicating missing pixels for later completion.
	- Not compatible with Weighted demosaicing and Weighted interpolation modes.
	
- New interpolation mode: '4DKernelReg' (M. LE PENDU)
	- Enable/Disable with DecodeOptions.ResampMethod='4DKernelReg' (default='fast').
	- Implemented based on the description in [7].
	- 4D kernel regression interpolation method (Performs demosaicing and resampling interpolations simultaneously).
	- Not compatible with Weighted demosaicing and Weighted interpolation modes.
	

## Other tools previously added in the version v0.4-CLIM:

- Weighted demosaicing (P. DAVID)
	- Enable/Disable with DecodeOptions.WeightedDemosaic (default=false).
	- Method detailed in [3].

- Weighted interpolation for rotation (P. DAVID)
	- Enable/Disable with DecodeOptions.WeightedInterp (default=true).
	- Method detailed in [3].
	- Not implemented for barycentric interpolation.

- Barycentric interpolation (method from KAIST toolbox [4], added by M. LE PENDU)
	- Enable/Disable with DecodeOptions.ResampMethod='barycentric' (default='fast').
	- Slower demultiplexing but generates higher resolution images.
	- Does not include generation of Weights (i.e. demultiplexed white image).

- Hot pixel correction (method from KAIST toolbox [4], added by M. LE PENDU)
	- Enable/Disable with DecodeOptions.HotPixelCorrect (default=true).
	- Note: Only works with F01 cameras (black images not available in ILLUM calibration data). For ILLUM data external hot pixel correction is necessary.

- Automatic White balancing (R. DAUDT)
	- Enable/Disable with DecodeOptions.DoAWB (default=false) and DecodeOptions.OptionalTasks='ColourCorrect' (default={}).
	- Method detailed in [5]. Original implementation from https://web.stanford.edu/~sujason/ColorBalancing/Code/robustAWB.m
	- Adapted for Light fields : White balance parameters determined for the central view and then applied to all the views.

- Alternative calibration using existing mlaCalib information from LYTRO metadata (S. HEMAMI, M. LE PENDU)
	- Called automatically by LFUtilProcessWhiteImages.m during calibration phase in addition to the calibration from the original toolbox.
	- To use mlaCalibration instead of the toolbox calibration for the decoding : set the decode option WhiteProcDataFnameExtension to '.grid2.json' when calling LFLytroDecodeImage (or LFUtilDecodeLytroFolder). (by default, '.grid.json' is set to use the toolbox calibration).
	- Warning : there may be errors in the way to derive the lenslet grid models (stored in .grid2.json files) from the mla calib data.


## Usage

Usage is similar to the Light field toolbox v0.4 with additional/modified DecodeOptions:

- .ResampMethod = 'fast' (default) / 'triangulation' / 'barycentric' / 'none' / '4DKernelReg'
	- option defined in LFDecodeLensletImageSimple.m

- .HotPixelCorrect = true (default) / false
	- option defined in LFLytroDecodeImage.m

- .WeightedDemosaic = true / false (default)
	- option defined in LFDecodeLensletImageSimple.m

- .WeightedInterp = true (default) / false
	- option defined in LFDecodeLensletImageSimple.m

- .CorrectSaturated = true (default) / false
	- option defined in LFUtilDecodeLytroFolder.m

- .noClip = true / false (default)
	- option defined in LFUtilDecodeLytroFolder.m

- .EarlyWhiteBalance = true (default) / false
	- option defined in LFUtilDecodeLytroFolder.m

- .DoAWB = true (use Automatic White Balancing) / false (default).
	- option defined in LFUtilDecodeLytroFolder.m
	- only applies if DecodeOptions.OptionalTasks='ColourCorrect'.

- .WhiteProcDataFnameExtension = '.grid.json' (default, as in LF toolbox v0.4) / '.grid2.json' (use lensletGridModels computed from mlaCalibration)
	- option defined in LFLytroDecodeImage.m

(The results in [6] were generated using DecodeOptions.OptionalTasks='ColourCorrect'. For the other options, default values were used)


## References

[1] P. Matysiak, M. Grogan, M. Le Pendu, M. Alain and A. Smolic, ["A Pipeline for Lenslet Light Field Quality Enhancement"](https://v-sense.scss.tcd.ie/research/light-fields/high-quality-light-field-extraction/), International Conference on Image Processing (ICIP) 2018.

[2] D. G. Dansereau, O. Pizarro, and S. B. Williams, "Decoding, calibration and rectification for lenselet-based plenoptic cameras," in Computer Vision and Pattern Recognition (CVPR), 2013.

[3] P. David, M. Le Pendu and C. Guillemot, ["White Lenslet Image Guided Demosaicing for Plenoptic Cameras"](https://www.irisa.fr/temics/demos/lightField/Demosaicing/LensletDemosaicing.html),  IEEE International Workshop on Multimedia Signal Processing (MMSP) 2017.

[4] D. Cho, M. Lee, S. Kim and Y.-W. Tai, "Modeling the Calibration Pipeline of the Lytro Camera for High Quality Light-Field Image Reconstruction", IEEE International Conference on Computer Vision (ICCV) 2013.

[5] J.-Y. Huo, Y.-L. Chang, J. Wang and X.-X. Wei "Robust Automatic White Balance Algorithm using Gray Color Points in Images", IEEE Transactions on Consumer Electronics, 2006.

[6] P. Matysiak, M. Grogan, M. Le Pendu, M. Alain, E. Zerman and A. Smolic, ["High Quality Light Field Extraction and Post-Processing for Raw Plenoptic Data"](https://v-sense.scss.tcd.ie/research/light-fields/high-quality-light-field-extraction/), Transactions on Image Processing 2020.

[7] S. Xu, Z.-L. Zhou, and N. Devaney, "Multi-view Image Restoration from Plenoptic Raw Images", Asian Conference on Computer Vision (ACCV) Workshops, 2014.