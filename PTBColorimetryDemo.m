% Colorimetry Demo

% Gives a demo of some of the datasets and computations available within
% PTB to do colorimetry.

% Additional goal: Provide a standard MATLAB implementation of CIE
% colorimetry to make sure that when researchers say they're computing X
% there's a slightly better chance it's comparable across studies.

% See also:
% - RenderDemo.m
% - DKLDemo.m
% - PsychColorimetricData/PsychColorimetricMatFiles/Contents.m

% Also: (Recommendation from DB, not gone through yet)

% PsychDemos/CalDemo
% PsychDemos/IsomerizationsInDishDemo
% PsychDemos/IsomerizationsInEyeDemo
% PsychDemos/NomogramDemo
% PsychTests/CIEConeFundamentalsTest
% PsychTests/CIEXYZPhysTest
% PsychTests/CIEConeFundamentalsFieldSizeTest
% PsychTests/OSAUCSTest
% PsychTests/TrolandTest
% PsychColorimetric/PsychMunsell/MunsellTest
% PsychRadiometric/RadiometricConversionsTest
% https://github.com/BrainardLab/TeachingCode :
% - MatlabTutorials/calibrationTutorial
% - MatlabTutorials/colorSpaceExample
% - MatlabTutorials/colorTransformAssignment
% - MatlabTutorials/silentSubstitutionTutorial
% Various in isetbio

% TO DO:
% - Improve intro
% - Try to include something that requires luminance calcs
% - u'v' space
% - MB space
% - Colour constancy (?)
% - Reference CIE docs
% - Refer to Westland toolbox

clc, clear, close all

%% Define and load: illuminants, surfaces and observer

% Illuminant - aka Spectral Power Distribution (SPD)
% In this example we use D65 which is considered a standard daylight 
% (https://en.wikipedia.org/wiki/Illuminant_D65)
% For explanation of PTB colorimetric data standards see PsychColorimetricData/PsychColorimetricMatFiles/Contents.m
load spd_D65

% Surface - aka Surface Reflectance Function (SRF)
% Surfaces of:
% Vrhel, M.J., Gershon, R. and Iwan, L.S., 1994. Measurement and Analysis of Object Reflectance Spectra. Color Research & Application, 19(1), pp.4–9.
load sur_vrhel

% Observer - aka Spectral Sensitivity Function (SSF)
% CIE 1931 observer:
% https://en.wikipedia.org/wiki/CIE_1931_color_space
load T_xyz1931

% These are just examples of many options, see PsychColorimetricData/PsychColorimetricMatFiles/Contents.m

%% Interpolate values to match
% Both S_D65 and S_xyz1931 have range/intervals of [380,5,81], so we'll
% choose to modify sur_vrhel (and S_vrhel)

sur_vrhel = SplineSrf(S_vrhel,sur_vrhel,S_D65,1);
S_vrhel = S_D65;

figure, plot(SToWls(S_vrhel), sur_vrhel)
xlabel('Wavelength (nm)')
ylabel('Reflectance')
% Note how the flag at the end of the SplineSrf call has requested that
% extrapolation be done by extending with final values.

%% Compute colorimetry

colourSignals = sur_vrhel.*spd_D65;
XYZ = T_xyz1931*colourSignals;
xy = [XYZ(1,:)./sum(XYZ);XYZ(2,:)./sum(XYZ)]; % Or use XYZToxyY.m

%% Plot colorimetry

spectralLocus = [T_xyz1931(1,:)./sum(T_xyz1931);T_xyz1931(2,:)./sum(T_xyz1931)];
sRGBSpectralLocus = XYZToSRGBPrimary(T_xyz1931); %These values go considerably out of gamut, but we only want rough values to orient ourselves

figure, hold on, 
axis equal
axis([0 1 0 1])
scatter(spectralLocus(1,:),spectralLocus(2,:),[],sRGBSpectralLocus','filled')
scatter(xy(1,:),xy(2,:),'kv')
xlabel('x')
ylabel('y')

%% - %%

% That's the basic stuff.
% Now let's break out some of the fancier stuff...

%% Calculate CIELAB values

whiteXYZ = T_xyz1931*spd_D65;
Lab = XYZToLab(XYZ,whiteXYZ);

figure, 
scatter3(Lab(2,:),Lab(3,:),Lab(1,:)) %To plot L on the Z axis

%% Calculate the Correlated Colour Temperature of an illuminant
SPDToCCT(spd_D65,S_D65)
% Why 6504K? See https://en.wikipedia.org/wiki/Illuminant_D65#Why_6504_K?




