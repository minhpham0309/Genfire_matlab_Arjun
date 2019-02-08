%% computing voronoi area %%

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




function obj = voronoi_area(obj)

projections = obj.InputProjections;
angles = obj.InputAngles;
particleWindowSize = obj.particleWindowSize(2);
interpolationCutoffDistance = obj.interpolationCutoffDistance;
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
    kMeasured = zeros(floor(particleWindowSize*oversamplingRatio),floor(particleWindowSize*oversamplingRatio),size(projections,3),'double');
else
    kMeasured = zeros(1 + floor((particleWindowSize-1)*oversamplingRatio),1 + floor((particleWindowSize-1)*oversamplingRatio),size(projections,3),'double');
end

tic %start clock

%get the dimension (assumed square and even) and setup the center and radius of the array size
dim1 = size(kMeasured,1);
nc = double(round((dim1+1)/2));%center pixel
n2 = single(nc-1);%radius of array

%setup the coordinates of the reciprocal slice to determine its 3D coordinates
[ky, kx] = meshgrid((1:1:dim1)-nc,(1:1:dim1)-nc);
ky = single(ky);kx = single(kx);
%Q = sqrt(ky.^2+kx.^2)./n2;
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
    for projNum = 1:size(projections,3);
        kMeasured(:,:,projNum) = my_fft(My_paddzero(projections(:,:,projNum),[size(kMeasured,1) size(kMeasured,2)]));
        %kMeasured2(:,:,projNum) = fftshift(fftn(My_paddzero(projections(:,:,projNum),[size(kMeasured,1) size(kMeasured,2)])));
    end
    %{
    figure; img(projections);
    figure; img(kMeasured);
    figure; img(kMeasured2);
    drawnow();
    %}
end
clear projections

%initialize arrays to contain coordinates
measuredX = zeros(size(kMeasured,2)*size(kMeasured,1),size(kMeasured,3),'double');
measuredY = zeros(size(kMeasured,2)*size(kMeasured,1),size(kMeasured,3),'double');
measuredZ = zeros(size(kMeasured,2)*size(kMeasured,1),size(kMeasured,3),'double');


numProj = size(kMeasured,3);
for projNum = 1:numProj;
    phi = angles(projNum,1);
    theta = angles(projNum,2);
    psi = angles(projNum,3);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  GENFIRE/RELION/XMIPP/FREALIGN/EMAN Euler angle convention:
    % %
    R = [ cosd(psi)*cosd(theta)*cosd(phi)-sind(psi)*sind(phi) ,cosd(psi)*cosd(theta)*sind(phi)+sind(psi)*cosd(phi)   ,    -cosd(psi)*sind(theta);
        -sind(psi)*cosd(theta)*cosd(phi)-cosd(psi)*sind(phi), -sind(psi)*cosd(theta)*sind(phi)+cosd(psi)*cosd(phi) ,   sind(psi)*sind(theta)  ;
        sind(theta)*cosd(phi)                               , sind(theta)*sind(phi)                                ,              cosd(theta)];
    
    rotkCoords = R'*[kx;ky;kz];%rotate coordinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    if projNum==1
        theta
        R
    end
    %}
    measuredX(:,projNum) = rotkCoords(1,:);%rotated X
    measuredY(:,projNum) = rotkCoords(2,:);%rotated Y
    measuredZ(:,projNum) = rotkCoords(3,:);%rotated Z
end



%reshape to simplify
%measuredX = reshape(measuredX,1,size(kMeasured,2)*size(kMeasured,1)*size(kMeasured,3));
%measuredY = reshape(measuredY,1,size(kMeasured,2)*size(kMeasured,1)*size(kMeasured,3));
%measuredZ = reshape(measuredZ,1,size(kMeasured,2)*size(kMeasured,1)*size(kMeasured,3));
%kMeasured = reshape(kMeasured,1,size(kMeasured,2)*size(kMeasured,1)*size(kMeasured,3));
badInd = find(kMeasured==-999);%delete values that are flagged as bad
measuredX(badInd) = [];
measuredY(badInd) = [];
measuredZ(badInd) = [];
kMeasured(badInd) = [];


if all(angles(projNum,1)==0) && all(angles(projNum,3)==0)
    fprintf('2D voronoi diagram\n');
    lower = -nc+1; upper=dim1-nc;
    [X,Z] = meshgrid((1:dim1)-nc,(1:dim1)-nc);
    [X_bdr1,Z_bdr1] = meshgrid([lower-1,upper+1]  , lower-1:upper+1);
    [X_bdr2,Z_bdr2] = meshgrid(lower:upper, [lower-1,upper+1]);
    X_bdr = [X_bdr1(:);X_bdr2(:)];
    Z_bdr = [Z_bdr1(:);Z_bdr2(:)];
    
    %[min(measuredX(:)),max(measuredX(:))]
    %[min(X(:)),max(X(:))]
    
    %[length(unique(measuredX + 1j*measuredZ)),size(measuredX)]
    %[~, measured_index,~] = unique(measuredX(:) + 1j*measuredZ(:));
    measuredX = reshape(measuredX,[dim1,dim1,numProj]);
    measuredX = permute(measuredX,[1,3,2]);
    measuredX = reshape(measuredX,[dim1*numProj,dim1]);
    
    measuredY = reshape(measuredY,[dim1,dim1,numProj]);
    measuredY = permute(measuredY,[1,3,2]);
    measuredY = reshape(measuredY,[dim1*numProj,dim1]);
    
    measuredZ = reshape(measuredZ,[dim1,dim1,numProj]);
    measuredZ = permute(measuredZ,[1,3,2]);
    measuredZ = reshape(measuredZ,[dim1*numProj,dim1]);
    
    kMeasured = reshape(kMeasured,[dim1,dim1,numProj]);
    kMeasured = permute(kMeasured,[1,3,2]);
    kMeasured = reshape(kMeasured,[dim1*numProj,dim1]);    
    
    
    measuredX_s = measuredX(:,1);
    measuredZ_s = measuredZ(:,1);
    
    % check unique: need unique_ind for define kMeasured
    [unique_val,unique_ind,~]=unique(measuredX_s+1j*measuredZ_s);
    fprintf('number of unique values in sampling points\n');
    [length(unique_val),length(measuredX_s)]
    
    measuredX_s = measuredX_s(unique_ind);
    measuredZ_s = measuredZ_s(unique_ind);
    
    figure; scatter(measuredX_s,measuredZ_s,10,'filled'); drawnow
    
    
    % check whether sampling points overlap with rectangle grids  
    % mem_overlap = ismember(measuredX_s,X) & ismember(measuredZ_s,Z);
    [is_overlap,~] = ismember(measuredX_s + 1j*measuredZ_s,X+1j*Z);
    fprintf('number of overlap with grid = %d\n',sum(is_overlap(:)));
    measuredX_r = measuredX_s(~is_overlap);
    measuredZ_r = measuredZ_s(~is_overlap);
    
    %index_overlap = index_overlap(index_overlap>0);
    
    fprintf('overlap member\n');
    [measuredX_s(is_overlap),measuredZ_s(is_overlap)]
    
    size(measuredX)
    
    measuredX_o = measuredX(unique_ind(is_overlap),:); measuredX_o = measuredX_o(:);
    measuredY_o = measuredY(unique_ind(is_overlap),:); measuredY_o = measuredY_o(:);
    measuredZ_o = measuredZ(unique_ind(is_overlap),:); measuredZ_o = measuredZ_o(:);
    kMeasured_o = kMeasured(unique_ind(is_overlap),:); kMeasured_o = kMeasured_o(:);
    
    measuredX = measuredX(unique_ind(~is_overlap),:);
    measuredY = measuredY(unique_ind(~is_overlap),:);
    measuredZ = measuredZ(unique_ind(~is_overlap),:);    
    kMeasured = kMeasured(unique_ind(~is_overlap),:);
    
    size(measuredX)
    
    %figure; scatter3(measuredX(:,1),measuredZ(:,1),abs(kMeasured(:,1))); drawnow
    %figure; scatter3(measuredX(:,2),measuredZ(:,2),abs(kMeasured(:,2))); drawnow
    
    %[measuredX(:,1)+1j*measuredZ(:,1),measuredX(:,2)+1j*measuredZ(:,2),measuredY(:,1),measuredY(:,2)]

    measuredX = measuredX(:);
    measuredY = measuredY(:);
    measuredZ = measuredZ(:);
    kMeasured = kMeasured(:);
        
    
    %measuredX_ovl = measuredX_s(is_overlap);
    %measuredZ_ovl = measuredZ_s(is_overlap);
    
    %fprintf('overlap member is\n');
    %[measuredX_s(is_overlap),measuredZ_s(is_overlap)]
    
    %size(measuredZ_r)
    mem_overlap2 = ismember(measuredX_r,X_bdr) & ismember(measuredZ_r,Z_bdr);
    [measuredX_r(mem_overlap2),measuredZ_r(mem_overlap2)]
    
    data = [[X(:);measuredX_r;X_bdr],[Z(:);measuredZ_r;Z_bdr]];
    
    %tic; TRIA = delaunay(data);toc
    %fprintf('done delaunay diagram\n');
    
    %return
    tic;[v,c] = voronoin(data);toc
    fprintf('done voronoi diagram\n');
    len = size(c,1);
    index = (zeros(len,1))>0;
    coef_area = zeros(len,1);
    tic
    for i=1:len
        index(i) =  all(c{i}~=1);
        if index(i)
            coef_area(i) = polyarea(v(c{i},1),v(c{i},2));
        end
    end
    toc;
    %[len,sum(index(:))]
    coef_area = coef_area(index);
    %coef_area = cellfun(@(c) polyarea(v(c,1),v(c,2)), (c(index)),'uniformoutput',true);
    fprintf('done computing area\n');
    [sum(isinf(coef_area(:))),sum(isnan(coef_area(:)))]
    
    % coefficient on grid
    coef_grid = coef_area(1:dim1^2);
    coef_grid = reshape(coef_grid,[dim1,dim1])';
    coef_grid = coef_grid(:);
    

    figure; img(reshape(coef_grid,dim1,dim1)); drawnow();

    %coef_grid2d = reshape(coef_grid,dim1,dim1);save('coef_grid.mat','coef_grid2d');
    
    coef_grid = reshape(coef_grid,[dim1,dim1]);
    coef_grid = repmat(coef_grid,[1,1,dim1]);
    coef_grid = permute(coef_grid,[1,3,2]);
    coef_grid = coef_grid(:);
    fprintf('number of zero in coef_grid %d\n',length(find(coef_grid==0)));
    
    
    % coefficient off grid
    %measuredX_f = repmat(measuredX_r,[dim1,1]);
    %measuredZ_f = repmat(measuredZ_r,[dim1,1]);
    %measuredY_f = repmat((1:dim1)-nc,[dim1*projNum,1]);
    %measuredY_f = measuredY_f(:);
    coef_spl = coef_area(dim1^2+1:end);    
    %save('coef_spl.mat','coef_spl');
    coef_spl = repmat(coef_spl,[dim1,1]);
    coef_spl = coef_spl(:);
    fprintf('number of zero in coef_spl %d\n',length(find(coef_spl==0)));
    %[length(coef_spl),length(measuredX_f)]
    %[max(measuredX_f), max(measuredY_f), max(measuredZ_f)]
    %[min(measuredX_f), min(measuredY_f), min(measuredZ_f)]
    
    % grid switching XX and YY
    [YY,XX,ZZ] = meshgrid((1:dim1)-nc,(1:dim1)-nc,(1:dim1)-nc);
    XX = XX(:);
    YY = YY(:);
    ZZ = ZZ(:);
    %index_overlap = ismember(XX+1j*ZZ, measuredX_ovl+1j*measuredZ_ovl);
    %coef_grid(index_overlap)=0;
    
    index_o = zeros(size(measuredX_o));
    for i=1:length(measuredX_o)
        index_o(i) = find(XX==measuredX_o(i) & YY==measuredY_o(i) & ZZ==measuredZ_o(i));
    end
    %index_c1 = (measuredX_o+nc-1)*dim1^2 + (measuredY_o+nc-1)*dim1 + measuredZ_o+nc;
    %index_c2 = (measuredX_o+nc-1)*dim1^2 + (measuredZ_o+nc-1)*dim1 + measuredY_o+nc;
    %index_c3 = (measuredY_o+nc-1)*dim1^2 + (measuredX_o+nc-1)*dim1 + measuredZ_o+nc;
    %index_c4 = (measuredY_o+nc-1)*dim1^2 + (measuredZ_o+nc-1)*dim1 + measuredX_o+nc;
    %index_c5 = (measuredZ_o+nc-1)*dim1^2 + (measuredX_o+nc-1)*dim1 + measuredY_o+nc;
    index_c6 = (measuredZ_o+nc-1)*dim1^2 + (measuredY_o+nc-1)*dim1 + measuredX_o+nc;
    %norm_index = norm(index_o-index_c1)
    %norm_index = norm(index_o-index_c2)
    %norm_index = norm(index_o-index_c3)
    %norm_index = norm(index_o-index_c4)
    %norm_index = norm(index_o-index_c5)
    norm_index = norm(index_o-index_c6)
    
    % input
    max_j = ceil(dim1/2);
    xj = [XX;measuredX]/max_j*pi;
    yj = [YY;measuredY]/max_j*pi;
    zj = [ZZ;measuredZ]/max_j*pi;
    fk = [randn(dim1^3,1);kMeasured];
    coef = [coef_grid;coef_spl];
    
    [max(abs(xj)),max(abs(yj)),max(abs(zj))]
    [max(kMeasured),min(kMeasured),sum(isnan(kMeasured(:)))]
    [max(coef),min(coef),sum(isnan(coef(:)))]
    
    %histogram(coef); drawnow();
    nj = length(xj);
    [nj,length(fk),length(coef)]
    
    ratio = length(fk)/dim1^3;
    paddedSupport = My_paddzero(obj.Support, [dim1,dim1,dim1]) ~=0 ;
    %size(paddedSupport)
    %figure;img(paddedSupport);drawnow();
    
    
    %return
    tic
    for i=1:5
        u = nufft3d1(nj,xj,yj,zj,fk.*coef,1,1e-12,dim1,dim1,dim1);
        fprintf('iteration %d\n',i)
        u = u*ratio;
        u = reshape(u,[dim1,dim1,dim1]);
        u = max(0,real(u)).* paddedSupport;
        fk(1:dim1^3) = fftshift(fftn(ifftshift(u)));
        fk(index_o) = kMeasured_o;
    end
    toc
    u = My_stripzero(u,[obj.Dim1,obj.Dim1,obj.Dim1]);
    figure;img(permute(u,[2,3,1])); 
    figure;img(permute(u,[1,3,2])); 
    figure;img(u); 
    
else %%%%%%%%%%%%%%%%%%  3D case   %%%%%%%%%%%%%%%%%%%
    lower = -nc; upper=dim1-nc+1;
    [X,Y,Z] = meshgrid((1:dim1)-nc,(1:dim1)-nc,(1:dim1)-nc);
    [X_bdr1,Y_bdr1,Z_bdr1] = meshgrid([lower,upper], lower:upper, lower:upper);
    [X_bdr2,Y_bdr2,Z_bdr2] = meshgrid(lower+1:upper-1, [lower,upper], lower:upper);
    [X_bdr3,Y_bdr3,Z_bdr3] = meshgrid(lower+1:upper-1, lower+1:upper-1, [lower,upper]);
    X_bdr = [X_bdr1(:);X_bdr2(:);X_bdr3(:)];
    Y_bdr = [Y_bdr1(:);Y_bdr2(:);Y_bdr3(:)];
    Z_bdr = [Z_bdr1(:);Z_bdr2(:);Z_bdr3(:)];
    
    
    is_overlap = ismember(measuredX,X) & ismember(measuredZ,Z) & ismember(measuredY,Y);
    measuredX_r = measuredX(~is_overlap);
    measuredY_r = measuredY(~is_overlap);
    measuredZ_r = measuredZ(~is_overlap);
    
    data = [[X(:);measuredX_r';X_bdr],[Y(:);measuredY_r';Y_bdr],[Z(:);measuredZ_r';Z_bdr]];
    size(data)
    
    tic; TRIA = delaunay(data);toc
    fprintf('done delaunay diagram\n')
    %return
    [v,c] = voronoin(data);
    fprintf('done voronoi diagram');
    len = size(c,1)
    index = (zeros(len,1))>0;
    coef_area = zeros(len,1);
    parfor i=1:len
        index(i) =  all(c{i}~=1);
        coef_area(i) = polyarea(v(c{i},1),v(c{i},2));
    end
    %coef_area = cellfun(@(c) polyarea(v(c,1),v(c,2)), (c(index)),'uniformoutput',true);
    fprintf('done computing area');
end

return

masterInd = [];%masterInd will be a large list of the grid indices
masterVals = [];%complex values to include in weighted averaging for those grid points
masterDistances = [];%distance from measured value to grid point
% masterConfidenceWeights = [];

if allowMultipleGridMatches
    shiftMax = round(interpolationCutoffDistance);
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

for Yshift = -shiftMax:shiftMax
    for Xshift = -shiftMax:shiftMax
        for Zshift = -shiftMax:shiftMax
            tmpX = (round(measuredX)+Xshift); % apply shift
            tmpY = (round(measuredY)+Yshift);
            tmpZ = (round(measuredZ)+Zshift);
            tmpVals = kMeasured;
            distances = sqrt(abs(measuredX-tmpX).^2+abs(measuredY-tmpY).^2+abs(measuredZ-tmpZ).^2); %compute distance to nearest voxel
            tmpY = tmpY+nc; %shift origin
            tmpZ = tmpZ+nc;
            tmpX = tmpX+nc;
            goodInd = (~(tmpX>dim1|tmpX<1|tmpY>dim1|tmpY<1|tmpZ>dim1|tmpZ<1)) & distances<=interpolationCutoffDistance;%find candidate values
            masterInd = [masterInd sub2ind([dim1 dim1 dim1],tmpX(goodInd),tmpY(goodInd),tmpZ(goodInd))]; %append values to lists
            masterVals = [masterVals tmpVals(goodInd)];
            masterDistances = [masterDistances distances(goodInd)];
            
        end
    end
end

%clear measuredX
%clear measuredY
%clear measuredZ
clear confidenceWeights

% Now that we have a list of the complex values to grid, their coordinates,
% and their distances from the nearest voxel, we want to reorganize the
% data so that all values matched to a given voxel are in the same place,
% so that the weighted sum can be computed. The number of values matched to
% each voxel can vary, and although one could use cell arrays for this
% purpose, they are quite slow. Instead, one can simply sort the indices,
% and then find the unique values by looking at the difference in
% consecutive elements.

masterInd = masterInd';
index_exact = masterDistances<1e-9;
index_k = masterInd(index_exact);
k_exact = masterVals(index_exact);
%masterDistances = masterDistances + 1e-9;
masterDistances(~index_exact) = 1 ./ masterDistances(~index_exact);
%masterDistances = 1 ./ masterDistances;
%masterDistances(masterDistances>0) = (max(0,interpolationCutoffDistance-masterDistances(masterDistances>0))/interpolationCutoffDistance./masterDistances(masterDistances>0)).^2;
%masterDistances(isnan(masterDistances)) = 0;

%%%%%-nearest neighbor
%obj.measuredK = accumarray(masterInd,masterVals,[dim1^3 1],@min);

%obj.measuredK = accumarray(masterInd,masterVals.*masterDistances,[dim1^3 1],@sum_min3);
%sumWeights = accumarray(masterInd,masterDistances,[dim1^3 1]);
%
obj.measuredK = accumarray(masterInd,masterVals.*masterDistances,[dim1^3 1]);
sumWeights = accumarray(masterInd,masterDistances,[dim1^3 1]);
obj.measuredK(sumWeights>0) = obj.measuredK(sumWeights>0) ./ sumWeights(sumWeights>0);
%}
obj.measuredK(index_k) = k_exact;
%nnz(obj.measuredK)

%
nc=double(nc);
[XX,YY,ZZ] = meshgrid((1:dim1)-nc,(1:dim1)-nc,(1:dim1)-nc);
index_interp = sumWeights>0;
XX=XX(index_interp);YY=YY(index_interp);ZZ=ZZ(index_interp);
%tic;delaunayTriangulation(double(measuredX'),double(measuredY'),double(measuredZ'));toc
measuredK = griddata(double(measuredX),double(measuredY),double(measuredZ),double(kMeasured),YY,XX,ZZ);
measuredK(isnan(measuredK))=0;
obj.measuredK = zeros(dim1,dim1,dim1,'single');
%nnz(isnan(obj.measuredK))
obj.measuredK(index_interp)=measuredK;
obj.measuredK(index_k) = k_exact;
%obj.measuredK(sumWeights==0)=0;
%nnz(obj.measuredK)
%}
%%%%%% reshape and make Hermitian matrix
obj.measuredK = reshape(obj.measuredK,[dim1 dim1 dim1]);
obj.measuredK = hermitianSymmetrize(obj.measuredK);

obj.recIFFT = My_stripzero(real(my_ifft(obj.measuredK)), [obj.Dim1 obj.Dim2 obj.Dim1]);
timeTakenToFillInGrid = toc;
timeTakenToFillInGrid = round(10*timeTakenToFillInGrid)./10;
fprintf('GENFIRE: Fourier grid assembled in %.12g seconds.\n\n',timeTakenToFillInGrid);


end
