# Skin-Pigmentation-Calibration

MATLAB script repository for Cherenkov skin pigmentation calibration alogrothim 

# System Requirements
1) MATLAB - commerically available through MathWorks. Run and tested on Version 2023B. Currently oeprating on macOS Sonomoa 14.0. 
2) CDose - commercially avaialble through DoseOptics. Run and tested on Windows 10, Lenovo Thinkpad P51 Ensure that your PC meets the minimum specs (16GB ram, 500GB SSD, i7 processor, NVIDIA GTX 1050 Ti or better, Windows 10 Pro 64-bit, USB3). In order to acquire images, a BeamSite camera (DoseOptics), optial fiber and optical repeater will be required to connect to aquisition computer.

# User Guide
1) generate_SurfaceCT --> imports an anonymized patient CT file (folder containing DICOMs or .tif) and creates a surface-weighted, 2D projection in the view of the camera.
2) import_Plan_Cherenkov --> imports Cherenkov images and dose map images (from CDose) and scales each to the proper range. Then, individual patient data is dose-normalized by matching the dose per beam images to the Cherenkov per beam images
3) surfaceMaps --> generates 3D point cloud of patient CT in world coordinates relative to the camera position and finds the normal vectors from the patient surface that are parallel to the camera-to-isocenter vector
4) CT_correction --> generates patient-specific CT calibration factors based on their surface-weighted CT maps
5) process_allColor --> reads in all patient color images, demosaics, scales, and color corrects
6) bayer_filt --> called by process_allColor
7) Lcolor_mean --> take region of interest (ROI) in a color image to get the average R, G, and B values. The output is the average luminance value, which is a function of RGB,  L = 0.2126 * R + 0.7152 * G + 0.0722 * B
8) roi_mean --> mean value of an ROI on a 2D image
