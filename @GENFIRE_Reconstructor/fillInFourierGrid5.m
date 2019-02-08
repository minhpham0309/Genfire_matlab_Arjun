%%  fillInFourierGrid %%

%%inputs:
%%  projections - measured projections
%%  angles - Euler angles in the form 3xN_projections, where each projection has 3 angles in the form [phi;theta;psi]
%%  interpolationCutoffDistance - radius of sphere in which to include measured
%%  oversamplingRatio - oversampling ratio for the projections in each direction
%%      values when filling in grid. All points within this sphere will be weighted
%%      linearly by their inverse distance.
%%  interpolationCutoffDistance - radius of interpolation kernel
%%  doCTFcorrection - flag to correct for Contrast Transfer Function (CTF) in projections, requires CTFparameters
%%  CTFparameters - structure containing defocus values and defocus angle for each projection
%%  allowMultipleGridMatches - whether or not to allow each measured datapoint to be matched to multiple grid points

%%outputs:
%%  rec - inverse FFT of the assembled Fourier grid
%%  measuredK -assembled Fourier Grid

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.




function obj = fillInFourierGrid5(obj)

projections = obj.InputProjections;
angles = obj.InputAngles;
particleWindowSize = obj.particleWindowSize(2);
cutoff_dist = obj.interpolationCutoffDistance;
doCTFcorrection = obj.doCTFcorrection;
CTFparameters = obj.CTFparameters;
allowMultipleGridMatches = obj.allowMultipleGridMatches;
oversamplingRatio = obj.oversamplingRatio;
%,particleWindowSize,oversamplingRatio,interpolationCutoffDistance,doCTFcorrection, [], allowMultipleGridMatches

%create empty CTF parameters if not doing CTF correction
if ~doCTFcorrection
    CTFparameters = [];
end
if doCTFcorrection && nargin < 6
    error('GENFIRE: doCTFcorrection is turned on, but CTFparameters was not provided.\n\n')
end

%calculate padding parameters for the inputted window size
padding = round(particleWindowSize*(oversamplingRatio-1)/2);
centralPixel = size(projections,2)/2+1;
halfWindowSize = particleWindowSize/2;

%initialize array to hold measured data
if mod(particleWindowSize,2)==0
    kMeasured = zeros(floor(particleWindowSize*oversamplingRatio),floor(particleWindowSize*oversamplingRatio),size(projections,3),'single');
else
    kMeasured = zeros(1 + floor((particleWindowSize-1)*oversamplingRatio),1 + floor((particleWindowSize-1)*oversamplingRatio),size(projections,3),'single');
end

tic %start clock

%get the dimension (assumed square and even) and setup the center and radius of the array size
dim1 = size(kMeasured,1);
nc = single(round((dim1+1)/2));%center pixel
n2 = single(nc-1);%radius of array

%setup the coordinates of the reciprocal slice to determine its 3D coordinates
[ky, kx] = meshgrid((1:dim1)-nc,(1:dim1)-nc);ky = single(ky);kx = single(kx);
Q = sqrt(ky.^2+kx.^2)./n2;
kx = single(kx(:))'; ky = single(ky(:))'; %initialize coordinates of unrotate projection slice
kz = zeros(1,dim1*dim1,'single'); %0 degree rotation is a projection onto the X-Y plane, so all points have kz=0;

%check for the presence of some of the CTF correction options and set defaults if they are absent
if doCTFcorrection
    if isfield(CTFparameters,'CTFThrowOutThreshhold')
        CTFThrowOutThreshhold = CTFparameters(1).CTFThrowOutThreshhold; %value below which to not grid points that were suppressed by the CTF
    else
        CTFThrowOutThreshhold = 0.05;%default value
    end
    if isfield(CTFparameters,'ignore_first_peak')
        ignore_first_peak =  CTFparameters(1).ignore_first_peak;
    else
        ignore_first_peak = 0;
    end
    
    for projNum = 1:size(projections,3);
        %get Contrast Transfer Function (CTF)
        pjK = projections(:,:,projNum);
        centralPixelK = size(pjK,2)/2+1;
        
        %crop out the appropriate window
        pjK = pjK(centralPixelK-halfWindowSize:centralPixelK+halfWindowSize-1,centralPixelK-halfWindowSize:centralPixelK+halfWindowSize-1);%window projection
        
        pjK = my_fft(padarray(pjK,[padding padding 0]));%pad and take FFT
        
        %get the CTF
        [CTF, gamma] = ctf_correction(pjK,CTFparameters(projNum).defocusU,CTFparameters(projNum).defocusV,CTFparameters(projNum).defocusAngle,ignore_first_peak);%get CTF
        if CTFparameters(projNum).phaseFlip %this should always be on unless your projections have already been CTF corrected elsewhere
            pjK(CTF<0) = -1*pjK(CTF<0);%phase flip
        end
        
        if CTFparameters(projNum).correctAmplitudesWithWienerFilter
            
            %get dimensions of the CTF array
            dim1_2 = size(CTF,1);
            nc2 = single(round((dim1_2+1)/2));%center pixel
            n22 = single(nc2-1);%radius of array
            
            %reciprocal indices
            [ky2, kx2] = meshgrid(-n22:n22-1,-n22:n22-1);ky2 = single(ky2);kx2 = single(kx2);
            Q2 = sqrt(ky2.^2+kx2.^2)./n22;
            
            SSNR = ones(size(Q2));%initialize SSNR map
            %interpolate the SSNR array from the provided values of the SSNR per frequency shell
            SSNR(:) = interp1(linspace(0,1+1e-10,size(CTFparameters(projNum).SSNR,2)),CTFparameters(projNum).SSNR,Q2(:),'linear');%make weighting map from average FRC
            SSNR(isnan(SSNR)) = 0;
            wienerFilter = abs(CTF)./(abs(CTF).^2+(1./SSNR));%construct Wiener filter for CTF amplitude correction
            pjK = pjK.*wienerFilter;
        elseif CTFparameters(projNum).multiplyByCTFabs%multiplying by CTF boosts SNR and is most useful for datasets that are extremely noisy
            pjK = pjK.*abs(CTF);
        end
        
        if CTFThrowOutThreshhold>0 %recalculate CTF at new array size for throwing out values that were near CTF 0 crossover
            pjK(abs(CTF)<CTFThrowOutThreshhold & (gamma>(pi/2))) = -999;%flag values where CTF was near 0 to ignore for gridding, but ignore out to first peak
        end
        
        kMeasured(:,:,projNum) = pjK;
    end
    
    if CTFThrowOutThreshhold > 0     %flag values below the where the CTF was smaller than the CTFThrowOutThreshhold
        for projNum = 1:size(projections,3);
            pjK = projections(centralPixelK-halfWindowSize:centralPixelK+halfWindowSize-1,centralPixelK-halfWindowSize:centralPixelK+halfWindowSize-1,projNum);
            pjK = my_fft(padarray(pjK,[padding padding 0]));
            CTF = ctf_correction(pjK,CTFparameters(projNum).defocusU,CTFparameters(projNum).defocusV,CTFparameters(projNum).defocusAngle,ignore_first_peak);%get CTF
            pjK(abs(CTF)<CTFThrowOutThreshhold) = -999;%flag values where CTF was near 0 to ignore for gridding
            kMeasured(:,:,projNum) = pjK;
        end
    end
else
    %otherwise, add the projection to the stack of data with no further corrections
    %kMeasured2 = zeros(size(kMeasured),'single');
    [dimk1,dimk2,~] = size(kMeasured);
    nk = floor(dimk1/2)+1;
    for projNum = 1:size(projections,3);        
        prj_i = My_paddzero(projections(:,:,projNum),[size(kMeasured,1) size(kMeasured,2)]);
        kMeasured(:,:,projNum) = my_fft(prj_i);
        %kMeasured2(:,:,projNum) = fftshift(fftn(My_paddzero(projections(:,:,projNum),[size(kMeasured,1) size(kMeasured,2)])));
        %kMeasured_i = sum(bsxfun(@times, prj_i(:), exp(-1*1i*2*pi*(XX*index_r/dimk1+YY*index_r/dimk2))),1);
        %size(kMeasured_i)
        %kMeasured(:,:,projNum)=reshape(kMeasured_i,[dimk1,dimk2]);
    end
    fprintf('size of kMeasure\n');
    size(kMeasured)
    %{
    figure; img(projections);
    figure; img(kMeasured); 
    figure; img(kMeasured2); 
    drawnow();
    %}
end
%return
clear projections

%initialize arrays to contain coordinates
measuredX = zeros(size(kMeasured,2)*size(kMeasured,1),size(kMeasured,3),'single');
measuredY = zeros(size(kMeasured,2)*size(kMeasured,1),size(kMeasured,3),'single');
measuredZ = zeros(size(kMeasured,2)*size(kMeasured,1),size(kMeasured,3),'single');


for projNum = 1:size(kMeasured,3);
    phi = angles(projNum,1);
    theta = angles(projNum,2);
    psi = angles(projNum,3);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  GENFIRE/RELION/XMIPP/FREALIGN/EMAN Euler angle convention:
    % %
    R = [ cosd(psi)*cosd(theta)*cosd(phi)-sind(psi)*sind(phi) ,cosd(psi)*cosd(theta)*sind(phi)+sind(psi)*cosd(phi)   ,   -cosd(psi)*sind(theta);
        -sind(psi)*cosd(theta)*cosd(phi)-cosd(psi)*sind(phi), -sind(psi)*cosd(theta)*sind(phi)+cosd(psi)*cosd(phi) ,   sind(psi)*sind(theta)  ;
        sind(theta)*cosd(phi)                               , sind(theta)*sind(phi)                                ,              cosd(theta)];
    
    rotkCoords = R'*[kx;ky;kz];%rotate coordinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    measuredX(:,projNum) = rotkCoords(1,:);%rotated X
    measuredY(:,projNum) = rotkCoords(2,:);%rotated Y
    measuredZ(:,projNum) = rotkCoords(3,:);%rotated Z
end

%reshape to simplify
badInd = find(kMeasured==-999);%delete values that are flagged as bad

if allowMultipleGridMatches
    shiftMax = round(cutoff_dist);
else
    shiftMax = 0;
end
%The nearest grid point to a measured value can be found by rounding, but
%there can be more than one grid point within the cutoff sphere, so must
%search locally for other possibilities. However in practice I have found
%this search can slow the program down greatly, without significant change
%in final result. Even searching 1 voxel in either direction increases the
%number of calculations by 3^3 = 27; For this reason I have set shiftMax = 0 and
%just assign values to their closest voxel.



%%
%
sigma = obj.sigma;
s_size = 20;
rangeZ = (1-nk):s_size:(dimk1-nk); remainder = dimk1-nk - rangeZ(end);
if remainder>s_size/2, rangeZ = [rangeZ,dimk1-nk];
else rangeZ(end) = dimk1-nk;
end
lenZ = length(rangeZ);
measuredX = measuredX(:);
measuredY = measuredY(:);
measuredZ = measuredZ(:);
kMeasured = kMeasured(:);
[measuredXYZ,index_XYZ_unique,index_XYZ] = unique([measuredX,measuredY,measuredZ],'rows');
measuredX = measuredXYZ(:,1);
measuredY = measuredXYZ(:,2);
measuredZ = measuredXYZ(:,3);
kMeasured = accumarray(index_XYZ, kMeasured, [] ,@mean);
masterInd=[];
masterVals = [];
%{
if obj.DFT_doGPU
    measuredX = gpuArray(measuredX);
    measuredY = gpuArray(measuredY);
    measuredZ = gpuArray(measuredZ);
    kMeasured = gpuArray(kMeasured);
end
%}
for t=1:(lenZ-1)^3
    if mod(t,(lenZ-1)^2)==1, fprintf('\n   %d\n',ceil(t/(lenZ-1)^2));end
    %i = ceil(t/(lenX-1)^2);
    %t_temp = mod(t,(lenX-1)^2);
    %j = floor(t_temp/(lenX-1))+1;
    %k = mod(t_temp,lenX-1)+1;
    
    [i,j,k] = ind2sub([lenZ-1, lenZ-1, lenZ-1],t);
    %if t<(lenX-1)^2, disp([i,j,k]);end
    
    %index_ijk = find(measuredX>=rangeX(i) & measuredX<=rangeX(i+1) & measuredY>=rangeX(j) & measuredY<=rangeX(j+1) & measuredZ>=rangeX(k) & measuredZ<=rangeX(k+1));
    index_ijk = find((measuredX>=rangeZ(i) & measuredX<=rangeZ(i+1) & measuredZ>=rangeZ(j) & measuredZ<=rangeZ(j+1) & measuredY>=rangeZ(k) & measuredY<=rangeZ(k+1)));
    if mod(t,lenZ-1)==1,fprintf('%d.length = ',j);end
    if length(index_ijk)>0
        fprintf('%d, ',length(index_ijk));
        
        measuredX_i = measuredX(index_ijk);
        measuredY_i = measuredY(index_ijk);
        measuredZ_i = measuredZ(index_ijk);
        measuredk_i = kMeasured(index_ijk);
        
        rx = bsxfun(@minus, measuredX_i,measuredX_i');
        ry = bsxfun(@minus, measuredY_i,measuredY_i');
        rz = bsxfun(@minus, measuredZ_i,measuredZ_i');
        A = exp(-sigma^2*sqrt(rx.^2 + ry.^2 + rz.^2));        
        phi = A\measuredk_i;
        %A(A<exp(-2*sigma^2))=0;
        %A = sparse(double(A));
        %[phi,~] = pcg(A,measuredk_i,1e-12);
        
        tmpX_i = round(measuredX_i);
        tmpY_i = round(measuredY_i);
        tmpZ_i = round(measuredZ_i);
        dist = sqrt(abs(measuredX_i-tmpX_i).^2+abs(measuredY_i-tmpY_i).^2+abs(measuredZ_i-tmpZ_i).^2);
        goodInd = ~(tmpX_i>dim1-nc| tmpX_i<1-nc| tmpY_i>dim1-nc| tmpY_i<1-nc| tmpZ_i>dim1-nc| tmpZ_i<1-nc) & dist < cutoff_dist;
        tmpX_i = tmpX_i(goodInd);
        tmpY_i = tmpY_i(goodInd);
        tmpZ_i = tmpZ_i(goodInd);
        masterInd_i = sub2ind([dim1 dim1 dim1],tmpX_i+nc,tmpY_i+nc,tmpZ_i+nc);
        masterInd = [masterInd;gather(masterInd_i)];
        
        rx = bsxfun(@minus, tmpX_i, measuredX_i');
        ry = bsxfun(@minus, tmpY_i, measuredY_i');
        rz = bsxfun(@minus, tmpZ_i, measuredZ_i');
        A = exp(-sigma^2*sqrt(rx.^2 + ry.^2 + rz.^2));
        masterVal_i = A*phi;
        masterVals = [masterVals;gather(masterVal_i)];
        
        %clear A rx ry rz
    end
    if mod(t,lenZ-1)==0,fprintf('\n');end
end
masterVals(isnan(masterVals))=0;
obj.measuredK = accumarray(masterInd,masterVals,[dim1^3 1],@mean);
%}

%% reshape and make Hermitian matrix
obj.measuredK = reshape(obj.measuredK,[dim1 dim1 dim1]);
obj.measuredK = hermitianSymmetrize(obj.measuredK);

obj.recIFFT = My_stripzero(real(my_ifft(obj.measuredK)), [obj.Dim1 obj.Dim2 obj.Dim1]);
timeTakenToFillInGrid = toc;
timeTakenToFillInGrid = round(10*timeTakenToFillInGrid)./10;
fprintf('GENFIRE: Fourier grid assembled in %.12g seconds.\n\n',timeTakenToFillInGrid);


end

function s = sum_min3(x)
    [~,I] = sort(x);
    s = sum(x(I(1:min(3,length(x)))));
end