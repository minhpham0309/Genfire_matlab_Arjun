%% run refinement
prefixstr = 'NiPt_Mo_181118BA_t2';
addpath('functions')
addpath('nufftall-1.33');
%/u/project/miao/minhrose/tomography/GENFIRE_matlab/
% absolute path on hoffman2

%pj_filename              = 'data/vesicle_projections.mat';
%angle_filename           = 'data/vesicle_angles.mat';

pj_filename              = 'data/small3_projections.mat';
angle_filename           = 'data/small3_angles.mat';

%pj_filename              = 'data/NiPt_Mo_FS_171118_t2_Projs_hori_TrAnRlnd.mat';
%angle_filename           = 'data\ref_angles_NiPt_Mo_FS_171118_t2.mat';

results_filename         = 'data/result.mat';

%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                         %%
%%%                        Welcome to GENFIRE!                              %%
%%%           GENeralized Fourier Iterative REconstruction                  %%
%%%                                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Author: Alan (AJ) Pryor, Jr.
%%% email:  apryor6@gmail.com
%%% Jianwei (John) Miao Coherent Imaging Group
%%% University of California, Los Angeles
%%% Copyright (c) 2015. All Rights Reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~isdeployed
    %addpath ../src/
   % addpath ../data/
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          User Parameters                              %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PLOT_YN = 0;

GENFIRE = GENFIRE_Reconstructor();

%%% See the README for description of parameters

GENFIRE.filename_Projections = pj_filename;
GENFIRE.filename_Angles = angle_filename ;
GENFIRE.filename_Results = results_filename;

GENFIRE.numIterations = 50; 
GENFIRE.pixelSize = .5; 
GENFIRE.oversamplingRatio = 4;              % ratio = 4
GENFIRE.griddingMethod = 2;                 % griddingMethod=2 for DFT
GENFIRE.allowMultipleGridMatches = 1;
GENFIRE.constraintEnforcementMode = 3;      % no suppression
GENFIRE.interpolationCutoffDistance =.125;    % interpolation Cut-off distance
GENFIRE.DFT_doGPU = 0;
GENFIRE.DFT_CentroSymmetricity=0;
GENFIRE.sigma = 1.0;

GENFIRE.dt_type = 1;                        % dt
GENFIRE.ds = 1;                             % ds
GENFIRE.smooth = 0;

GENFIRE.constraintPositivity = 1;
GENFIRE.constraintSupport = 1;
GENFIRE.ComputeFourierShellCorrelation = 0; 
GENFIRE.numBins = 50;
GENFIRE.percentValuesForRfree = 0.05;
GENFIRE.numBinsRfree = 35;
GENFIRE.doCTFcorrection = 0;
GENFIRE.CTFThrowOutThreshhold = 0;
GENFIRE.calculate_Rfree = 0;



%GENFIRE.particleWindowSize = [];
%GENFIRE.phaseErrorSigmaTolerance = [];
%GENFIRE.vector1 = [0 0 1];
%GENFIRE.vector2 = [0 1 0];
GENFIRE.vector3 = [1 0 0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin GENFIRE

% Read data from files. The data can also be directly be fed into GENFIRE
% Class, not necessarily from file
GENFIRE = readFiles(GENFIRE);

% Check if the data is consistent, and prepare reconstruction variables
% based on given parameters
GENFIRE = CheckPrepareData(GENFIRE);

% Run GENFIRE gridding
GENFIRE = runGridding(GENFIRE); 
%GENFIRE = voronoi_area(GENFIRE); 

% Run FSC if flag is on
if GENFIRE.ComputeFourierShellCorrelation
  
    % Do not calculate R_free for FSC calculation
    calculate_Rfree_ori = GENFIRE.calculate_Rfree;
    GENFIRE.calculate_Rfree= 0;
    
    GENFIRE = runFSC(GENFIRE);
    
    % set the R_free flag back go original
    GENFIRE.calculate_Rfree=calculate_Rfree_ori;
end

%GENFIRE = ClearCalcVariables(GENFIRE);
    
%run reconstruction
%GENFIRE = reconstruct(GENFIRE);
GENFIRE = reconstruct_dr(GENFIRE);
reconstruction = GENFIRE.reconstruction;
%final_rec = GENFIRE.final_rec;
%%
[Rtotal,Rarr,simuProjs] = Tian_calc_Rfactor_realspace(obj.reconstruction,obj.InputProjections,obj.InputAngles);
%%
%reconstruction = finalvol;
min1 = min(abs(reconstruction(:))); max1 = max(abs(reconstruction(:)));
figure;img(permute(reconstruction,[2,3,1]),'caxis',[min1,max1]);
figure;img(permute(reconstruction,[1,3,2]),'caxis',[min1,max1]);
figure;img(reconstruction,'caxis',[min1,max1]);
%%
figure;img(sum(permute(reconstruction,[2,3,1]),3));
figure;img(sum(permute(reconstruction,[1,3,2]),3));
figure;img(sum(reconstruction,3));

%%
final_rec = GENFIRE.final_rec;
figure;img(permute(final_rec,[2,3,1]));
figure;img(permute(final_rec,[1,3,2]));
figure;img(final_rec);
%%
save('Ni_DR.mat','final_rec');
%%
load('data\Up_NiPt_1118t2_FS.mat')
%%
my_rec = tom_mrcread('vesicle.mrc');
%%
[SumX_1118t2_dr_dt01,SumY_1118t2_dr_dt01,SumZ_1118t2_dr_dt01,ESTvol_dr_dt01,peakX_1118t2_dr_dt01,peakY_1118t2_dr_dt01,peakZ_1118t2_dr_dt01,~] = Rot_FinalVol(abs(reconstruction),Up);

figure;img(ESTvol_dr_dt01,'caxis',[0,max(ESTvol_dr_dt01(:))])
figure; img(SumX_1118t2_dr_dt01,'caxis',[0,max(SumX_1118t2_dr_dt01(:))]);

%%
%reconstruction = obj.reconstruction;
%load('NiPt_Mo_FS_171118_t2_Projs_hori_TrAnRlnd.mat');
%load('ref_angles_NiPt_Mo_FS_171118_t2.mat');

%InputProjections = shiftData;
InputProjections = GENFIRE.InputProjections;
Angles = GENFIRE.InputAngles;
reconstruction = GENFIRE.reconstruction;
NumProjs = size(InputProjections,3);
for i=1:NumProjs
    pj = InputProjections(:,:,i);
    back_pj = calculate3Dprojection_RealSpaceinterp(reconstruction,Angles(i,1),Angles(i,2),Angles(i,3));
    Rarr(i) = sum( abs(abs(back_pj(:))-abs(pj(:)))) / sum(abs(pj(:)));
    %fprintf('Rfactor: %.5f \n',Rarr(i));
    %simuProjs = cat(3,simuProjs,back_pj);
end
errF = mean(Rarr)

            

