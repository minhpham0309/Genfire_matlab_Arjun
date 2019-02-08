function [SumX,SumY,SumZ,RotatedVol,PeakPosX,PeakPosY,PeakPosZ,U] = Rot_FinalVol(final_Rec,U)
% inputDir = 'D:\Google Drive\Group\Projects\NiPtMo\WorkPerDay\20171118\GENFIRE_ForNewRefine_OriginMethod\results\';

% load([inputDir 'NiPt_reconstruction.mat']);

% savefileName = [inputDir 'NiPt_reconstruction_ROT.mat'];

%%
%%%Transfer the reconstruction to ESTvol  (you have to know what you import!)
% ESTvol = obj.reconstruction;

% ESTvol = GENFIRE_parameters.reconstruction;

ESTvol = final_Rec; %% for version 1

% ESTvol = NiPt_reconstruction;

a = 3.92;
% resolution =0.33;    %%%in anglstrom/pixel
% resolution =0.469;
resolution =0.495;

ESTvol_size1 = size(ESTvol,1);
ESTvol_size2 = size(ESTvol,2);
ESTvol_size3 = size(ESTvol,3);


% Center the arrays and make mesh grid
Yarr = My_find_symm_indarr(ESTvol_size1);
Xarr = My_find_symm_indarr(ESTvol_size2);
Zarr = My_find_symm_indarr(ESTvol_size3);

% [X1, Y1, Z1] = meshgrid(Xarr, Yarr, Zarr);

%% orinet
if nargin == 1
U   = OrientFCC(ESTvol,        resolution, a/(2),        .03,          5,           8,      1,      1,     7); % lattice constant has divided by 2.
% U = OrientFCC(reconstruction,res,        lattice_const,mag_tolerance,ang_tolerance,rankmax,userank,padnum,box3size)
end
%%

invU = inv(U);

% take only rotation part of the orientation matrix
[U1, S, U2] = svd(U);
Up = U1*U2';

% Set up meshgrid and rotate coordinates
[X1, Y1, Z1] = meshgrid(Xarr, Yarr, Zarr);

Xshape = size(X1);
Yshape = size(Y1);
Zshape = size(Z1);

posvec = [X1(:), Y1(:), Z1(:)]';
rot_posvec = Up*posvec;

rotX = reshape(rot_posvec(1,:),Xshape);
rotY = reshape(rot_posvec(2,:),Yshape);
rotZ = reshape(rot_posvec(3,:),Zshape);

clear posvec;
clear rot_posvec;

% Cut to make interpolation splitted, save memory
cut_size = 50;

cutnum = floor(size(rotX,3)/cut_size);
cutrem = mod(size(rotX,3),cut_size);

ROTvol = zeros(size(rotX));
for i=1:cutnum
    cutindar = ((i-1)*50+1):i*50;
    ROTtemp = interp3(X1,Y1,Z1,ESTvol,rotX(:,:,cutindar),rotY(:,:,cutindar),rotZ(:,:,cutindar)); %express the points in the rotated coordinates
    %ROTtemp = interp3(X,Y,Z,Support,rotX(:,:,cutindar),rotY(:,:,cutindar),rotZ(:,:,cutindar));
    ROTvol(:,:,cutindar) = ROTtemp;    
end

if cutrem~=0
    cutindar = (cutnum*50+1):size(rotX,3);
    ROTtemp = interp3(X1,Y1,Z1,ESTvol,rotX(:,:,cutindar),rotY(:,:,cutindar),rotZ(:,:,cutindar));
    %ROTtemp = interp3(X,Y,Z,Support,rotX(:,:,cutindar),rotY(:,:,cutindar),rotZ(:,:,cutindar));
    ROTvol(:,:,cutindar) = ROTtemp;
end

RotatedVol = ROTvol;
RotatedVol(isnan(RotatedVol)) = 0;

% save(savefileName, 'RotatedVol','U');
%%
% RotatedVol_NiPt_171118_tomo1_reBGsub_0107Rot = rot90(rot90(RotatedVol_NiPt_171118_tomo1_reBGsub_0107));

[SumX,PeakPosX]=sumSlicesPeak(RotatedVol,1);
[SumY,PeakPosY]=sumSlicesPeak(RotatedVol,2);
[SumZ,PeakPosZ]=sumSlicesPeak(RotatedVol,3);
%%
% save('sum_NiPt_20171118_tomo1_1226_62pj_NewComRot.mat', 'RotatedVol_NiPt_20171118_tomo1_1226_62pj_NewComRot','U','sumX_NiPt_20171118_tomo1_1226_62pj_NewComRot','sumY_NiPt_20171118_tomo1_1226_62pj_NewComRot','sumZ_NiPt_20171118_tomo1_1226_62pj_NewComRot');

end