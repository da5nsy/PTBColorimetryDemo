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
% For explanation of PTB colorimetric data format standards see PsychColorimetricData/PsychColorimetricMatFiles/Contents.m
load spd_D65

% Surface - aka Surface Reflectance Function (SRF)
% Surfaces of:
% Vrhel, M.J., Gershon, R. and Iwan, L.S., 1994. Measurement and Analysis of Object Reflectance Spectra. Color Research & Application, 19(1), pp.4–9.
load sur_vrhel

% Observer - aka Colour Matching Function (CMF) (aka Spectral Sensitivity Function, when it is biologically plausible)
% CIE 1931 observer:
% https://en.wikipedia.org/wiki/CIE_1931_color_space
load T_xyz1931

% These are just examples of many options, see PsychColorimetricData/PsychColorimetricMatFiles/Contents.m

%% Plot data

figure,
plot(SToWls(S_D65),spd_D65)

figure,
plot(SToWls(S_vrhel),sur_vrhel)

figure,
plot(SToWls(S_xyz1931),T_xyz1931)

%% Interpolate values to match
% Both S_D65 and S_xyz1931 have range/intervals of [380,5,81], so we'll
% choose to downsample/interpolate sur_vrhel (and S_vrhel)

sur_vrhel = SplineSrf(S_vrhel,sur_vrhel,S_D65,1);
S_vrhel = S_D65;

figure, plot(SToWls(S_vrhel), sur_vrhel)
xlabel('Wavelength (nm)')
ylabel('Reflectance')
% (Note how the flag at the end of the SplineSrf call has requested that
% extrapolation be done by extending with final values.)

%% Compute colorimetry

colourSignals = sur_vrhel.*spd_D65; % This is nominally the light actually reaching the eye

XYZ = T_xyz1931*colourSignals; % These are 'tristimulus values'

xy = [XYZ(1,:)./sum(XYZ);XYZ(2,:)./sum(XYZ)]; % 'Chromaticity co-ordinates' (Alternative: use XYZToxyY.m)

%% Plot colorimetry

spectralLocus_xy = [T_xyz1931(1,:)./sum(T_xyz1931);T_xyz1931(2,:)./sum(T_xyz1931)];
sRGBSpectralLocus = XYZToSRGBPrimary(T_xyz1931); %These values go considerably out of gamut, but we only want rough values to orient ourselves

figure, hold on, 
axis equal
axis([0 1 0 1])
scatter(spectralLocus_xy(1,:),spectralLocus_xy(2,:),[],sRGBSpectralLocus','filled')
scatter(xy(1,:),xy(2,:),'kv')
xlabel('x')
ylabel('y')

%% Calculate u'v' chromaticity

% More perceptually uniform than xy

upvp = xyTouv(xy); %upvp is short for u prime v prime, because there is another colourspace called uv which is subtly different *sigh*
sRGBSpectralLocus_upvp = xyTouv(spectralLocus_xy);

% manually from XYZ:
upvp_manFromXYZ = [4*XYZ(1,:)./(XYZ(1,:) + 15*XYZ(2,:) + 3*XYZ(3,:));...
    9*XYZ(2,:)./(XYZ(1,:) + 15*XYZ(2,:) + 3*XYZ(3,:))];
    
% manually from xy:
upvp_manFromxy = [4*xy(1,:)./(-2*xy(1,:) + 12*xy(2,:) + 3);...
    9*xy(2,:)./(-2*xy(1,:) + 12*xy(2,:) + 3)];

% check:
% find the maximum absolute difference between the above different methods
max(max(abs(upvp-upvp_manFromxy)))
max(max(abs(upvp-upvp_manFromXYZ)))
max(max(abs(upvp_manFromxy-upvp_manFromXYZ))) 
% should be tiny numbers (e.g. 2e-16), probably rounding errors

figure, hold on, 
axis equal
axis([0 1 0 1])
scatter(sRGBSpectralLocus_upvp(1,:),sRGBSpectralLocus_upvp(2,:),[],sRGBSpectralLocus','filled')
scatter(upvp(1,:),upvp(2,:),'kv')
xlabel('u''') %to get an apostrophe you throw extra apostrophes at it until it behaves
ylabel('v''')

%% CIELAB

whiteXYZ = T_xyz1931*spd_D65;
Lab = XYZToLab(XYZ,whiteXYZ);

figure, 
scatter3(Lab(2,:),Lab(3,:),Lab(1,:),'kv') 
xlabel('a*')
ylabel('b*')
zlabel('L*')

%% CIELUV

whiteXYZ = T_xyz1931*spd_D65;
Luv = XYZToLuv(XYZ,whiteXYZ);

figure, 
scatter3(Luv(2,:),Luv(3,:),Luv(1,:),'kv') 
xlabel('u*')
ylabel('y*')
zlabel('L*')

%% MB
load T_cones_ss2.mat
load T_CIE_Y2.mat
%plot(SToWls(S_cones_ss2),T_cones_ss2)

T_c = SplineCmf(S_cones_ss2,T_cones_ss2,S_xyz1931); %resampling so that I can use the old sRGBs that I already calculated, and keep the appearance comparable accross diagrams
T_C = SplineCmf(S_CIE_Y2,T_CIE_Y2,S_xyz1931);
spectralLocus_MB = LMSToMacBoyn(T_c,T_c,T_C);

LMS = T_c*colourSignals;
ls = LMSToMacBoyn(LMS,T_c,T_C);

figure, hold on
scatter(spectralLocus_MB(1,:),spectralLocus_MB(2,:),[],sRGBSpectralLocus','filled')
scatter(ls(1,:),ls(2,:),'kv')
xlabel('{\itl}_{MB}');
ylabel('{\its}_{MB}');


%% - %%

% That's the basic stuff.
% Now let's break out some of the fancier stuff...



%% Calculate the Correlated Colour Temperature of an illuminant
SPDToCCT(spd_D65,S_D65)
% Why 6504K? See https://en.wikipedia.org/wiki/Illuminant_D65#Why_6504_K?




