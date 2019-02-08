%%  GENFIRE_iterate %%

%%inputs:
%%  numIterations - number of iterations to run
%%  initialObject - inital guess of object
%%  support - region of 1's and 0's defining where the reconstruction can exist (voxels that are 0 in support will be set to 0 in reconstruction)
%%  measuredK - measured Fourier points
%%  constraintIndicators - a flag value that is used to determine when a particular Fourier point is enforced, i.e. by resolution.
%%      The values to enforce at any particular iteration are determined by constraintEnforcementDelayIndicators.
%%  constraintEnforcementDelayIndicators - vector of values indicating the Fourier enforcement cutoff. Fourier grid points with a constraintIndicators value greater than
%%      or equal to the current constraintEnforcementDelayIndicators value will be enforced. The values in constraintEnforcementDelayIndicators are spread evenly over the number of iterations.
%%  R_freeInd_complex - indices of 5% of points in the highest resolution shell of measuredK that are being withheld from
%%          reconstruction and compared to after iteration.
%%  R_freeVals_complex - corresponding complex values at the indices in R_freeInd_complex
%%  enforce_positivity - whether or not to enforce positivity constraint
%%  enforce_support - whether or not to enforce support constraint

%%outputs:
%%  rec - reconstruction after iteration
%%  errK - reciprocal space error
%%  Rfree_complex - value of complex Rfree in each resolution shell vs iteration

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015. All Rights Reserved.

function obj = reconstruct_dr(obj)
dt_type = obj.dt_type;
ds =    obj.ds;
obj.measuredK = ifftshift(obj.measuredK);
obj.measuredK_mask = ifftshift(obj.measuredK_mask);

numIterations = obj.numIterations;
enforce_positivity = obj.constraintPositivity;
enforce_support = obj.constraintSupport;

% paddedSupport = My_paddzero(obj.Support, [obj.n1_oversampled obj.n2_oversampled obj.n1_oversampled]);
paddedSupport = My_paddzero(obj.Support, size(obj.measuredK)) ~=0 ;

Q = make_Kspace_indices(paddedSupport);
Q=ifftshift(Q);
resRange = -0.05;%thickness of resolution ring to use for removal of datapoints for Rfree test

constraintIndicators = zeros(size(Q)); %% first reconstruction uses resolution extension/suppression, where lower resolution information is enforced
%%initially and the maximum enforced resolution increases. This is followed by resolution where suppression, which is the same process run backwards, so at the final iteration
%%only the lowest resolution is being enforced again.
constraintIndicators(obj.measuredK~=0 & obj.measuredK_mask) = 1-Q(obj.measuredK~=0 & obj.measuredK_mask);%make lower resolution have higher confidence. Higher confidence means enforced earlier in the reconstruction and for longer overall than low confidence

%remove datapoints for Rfree calculation

if obj.calculate_Rfree==1
    spatialFrequencyForRfree = linspace(0,1,obj.numBinsRfree+1);%compute spatial frequency bins
    Rfree_complex_bybin = zeros(obj.numBinsRfree,numIterations);
    for shellNum = 1:obj.numBinsRfree %loop over each frequency shell
        
        measuredPointInd_complex = find(obj.measuredK~=0&Q>=(spatialFrequencyForRfree(shellNum)+resRange)&Q<spatialFrequencyForRfree(shellNum+1)); %candidate values for Rfree_complex
        
        if ~isempty(obj.RandomSeed)
            rng(obj.RandomSeed);
        end
        P = randperm(numel(measuredPointInd_complex)); %shuffle values
        
        measuredPointInd_complex = measuredPointInd_complex(P); %apply shuffle
        cutoffInd_complex = floor(numel(measuredPointInd_complex).*obj.percentValuesForRfree); %take indices for 5% of measured data
        if cutoffInd_complex == 0 %make sure to include at least one point
            cutoffInd_complex = 1;
        end
        R_freeInd_complex{shellNum} = measuredPointInd_complex(1:cutoffInd_complex);%take complex value for 5% of measured data
        %R_freeInd_shifted_complex{shellNum} = My_iffshift3_ind(size(paddedSupport),R_freeInd_complex{shellNum});
        %now create a temporary set of constraints that have this 5% of
        %datapoints removed
        obj.measuredK_mask(R_freeInd_complex{shellNum}) = false;
        R_freeVals_complex{shellNum} = obj.measuredK(R_freeInd_complex{shellNum});
    end
end

%run the actual reconstruction
fprintf('GENFIRE Douglas Rachford: Reconstructing... \n\n');

tic;

if isempty(obj.initialObject)
    initialObject = zeros(size(paddedSupport),'single');
else
    initialObject = My_paddzero(obj.InitialObject, [obj.n1_oversampled obj.n2_oversampled obj.n1_oversampled]);
end


% numIterations,initialObject,support,measuredK,constraintIndicators,constraintEnforcementDelayIndicators,
%
% R_freeInd_complex,R_freeVals_complex,

bestErr = 1e30;%initialize best error

if obj.calculate_Rfree==1
    obj.Rfree_complex = zeros(1,numIterations,'single');%% initialize Rfree_complex curve , -1 is a flag that means undefined
end
obj.errK = ones(1,numIterations,'single');
obj.errF = ones(1,numIterations,'single');

%prefetch indices to use for error metric to avoid having to lookup each
%iteration
errInd = obj.measuredK~=0&obj.measuredK_mask;
%errInd_shifted = int32(My_iffshift3_ind(size(paddedSupport),errInd));


%determine how to spread the provided weighting cutoffs over the iterations
iterationNumsToChangeCutoff = round(linspace(1,numIterations,numel(obj.constraintEnforcementDelayIndicators)));

currentCutoffNum = 1;

% do ifftshift and fftshift only at before and after iteration
% and avoid using ifftshift and fftshift during the iteration to save time
paddedSupport = ifftshift(paddedSupport);
initialObject = ifftshift(initialObject);

% single precision
initialObject = single(initialObject);
constraintIndicators = single(constraintIndicators);
clear Q
u = initialObject;

Vol_ind = obj.Vol_ind;
Rarr = zeros(obj.NumProjs,1);
Angles = (obj.InputAngles);

smooth = obj.smooth;
if smooth
    [d1,d2,d3] = size(initialObject);
    d1 = single(d1); d2 = single(d2); d3 = single(d3);
    [XX,YY,ZZ] = meshgrid(1:d1,1:d2,1:d3);
    x_cen = floor(d1/2);
    y_cen = floor(d2/2);
    z_cen = floor(d3/2);
    kernel = (XX-x_cen).^2 + (YY-y_cen).^2 + (ZZ-z_cen).^2;
    sigma = 1/obj.smooth;
    kernel = exp(-kernel/sigma^2);
    clear XX YY ZZ
end

for iterationNum = 1:numIterations
    if iterationNum == iterationNumsToChangeCutoff(currentCutoffNum) && iterationNum<2
        currentCutoffNum = find(iterationNumsToChangeCutoff==iterationNum,1,'last');
        %constraintInd_complex = (constraintIndicators>obj.constraintEnforcementDelayIndicators(currentCutoffNum)) & obj.measuredK~=0 & obj.measuredK_mask;
        %constraintInd_complex = obj.measuredK~=0 & obj.measuredK_mask;
        %constraintInd_complex_shifted = int32(My_iffshift3_ind(size(paddedSupport),constraintInd_complex));
        %constraintInd_complex_shifted = ifftshift(constraintInd_complex);
        currentCutoffNum = currentCutoffNum+1;
        bestErr = 1e30;%reset best error
    end
    
    switch dt_type
        case 0
            dt = 0;
        case 1
            dt = 0.1;
        case 2
            dt = 0.1 + 0.3*(1-sqrt((numIterations-iterationNum)/numIterations));
    end
    %ds = 0.5 * sqrt((numIterations-iterationNum)/numIterations) ;
    
    
    k = fftn(initialObject);%take FFT of initial object
    %monitor error
    obj.errK(iterationNum) = sum(abs(abs(k(errInd))-abs(obj.measuredK(errInd))))./sum(abs(obj.measuredK(errInd)));
    
    if obj.calculate_Rfree==1 %if values have been withheld from measuredK for monitoring R_free, check them accordingly
        if ~isempty(R_freeInd_complex)
            %calculate Rfree in each resolution shell
            total_Rfree_error      = 0;
            total_Rfree_error_norm = 0;
            for shellNum = 1:numel(R_freeInd_complex)
                %tmpInd =R_freeInd_complex{shellNum};
                tmpInd_shifted =R_freeInd_complex{shellNum};
                tmpVals = R_freeVals_complex{shellNum};
                obj.Rfree_complex(shellNum,iterationNum) = sum(abs(k(tmpInd_shifted)-tmpVals))./sum(abs(tmpVals));

                Rfree_numerator                         = sum(abs(k(tmpInd_shifted)-tmpVals));
                Rfree_denominator                       = sum(abs(tmpVals));
                total_Rfree_error                       = total_Rfree_error + Rfree_numerator;
                total_Rfree_error_norm                  = total_Rfree_error_norm + Rfree_denominator;
                Rfree_complex_bybin(shellNum, iterationNum) = Rfree_numerator / Rfree_denominator;
                
            end
            Rfree_complex_total(iterationNum) = total_Rfree_error / total_Rfree_error_norm;
        end
        %figure(100);
        %subplot(1,3,1);plot(mean(obj.Rfree_complex,1));title('Mean R-free Value v.s. Iteration ')
        %subplot(1,2,1);plot(Rfree_complex_total);xlim([0,numIterations]); title('Mean R-free Value v.s. Iteration ')
        %subplot(1,2,2);plot(Rfree_complex_bybin(:,iterationNum));title('Final R-free Values v.s. Spatial Frequency');
        %drawnow();
    end
    
    % computing R_real
    %
    if mod(iterationNum,10)==0
        for i=1:obj.NumProjs
            pj = obj.InputProjections(:,:,i);
            model = fftshift(initialObject);
            model = model(Vol_ind(1,1):Vol_ind(1,2),Vol_ind(2,1):Vol_ind(2,2),Vol_ind(3,1):Vol_ind(3,2));
            back_pj = calculate3Dprojection_RealSpaceinterp(model,Angles(i,1),Angles(i,2),Angles(i,3));                        
            Rarr(i) = sum( abs(abs(back_pj(:))-abs(pj(:)))) / sum(abs(pj(:)));
            %fprintf('Rfactor: %.5f \n',Rarr(i));
            %simuProjs = cat(3,simuProjs,back_pj);
        end
        errF = mean(Rarr);
        fprintf('Iteration %d. Rfactor: %.6f \n',iterationNum, errF);
        obj.errF(iterationNum) = errF;
    end
    %}    
    
    if obj.errF(iterationNum)<bestErr %if current reconstruction has better error, update best error and best reconstruction
        bestErr = obj.errK(iterationNum);
        obj.reconstruction = fftshift(initialObject);
        obj.reconstruction = obj.reconstruction(obj.Vol_ind(1,1):obj.Vol_ind(1,2), obj.Vol_ind(2,1):obj.Vol_ind(2,2), obj.Vol_ind(3,1):obj.Vol_ind(3,2));
    end
    
    
    fprintf('Iteration %d. Error = %.6f\n',iterationNum, obj.errK(iterationNum));
    %enforce Fourier constraint
    %k(constraintInd_complex_shifted) = obj.measuredK(constraintInd_complex);
    %k(constraintInd_complex) = dt*k(constraintInd_complex) + (1-dt)*obj.measuredK(constraintInd_complex);
    k(errInd) = dt*k(errInd) + (1-dt)*obj.measuredK(errInd);
    u_K = ifftn(k);
    %u_K = real(ifftn(k));
    initialObject = (1+ds)*u_K - ds*u;
    
    %     initialObject = real(ifftn(k));%obtain next object with IFFT
    %obj_temp = initialObject;
    initialObject = real(initialObject);%obtain next object with IFFT
    initialObject = max(0,initialObject);
    initialObject = initialObject.*paddedSupport;
    if smooth 
        %if mod(iterationNum,20)==0, kernel = kernel.^1.1 ;end
        F_obj = fftn(initialObject) .* kernel; initialObject = real(ifftn(F_obj));
        %F_obj = fftn(fftshift(initialObject)) .* kernel; initialObject = real(ifftshift(ifftn(F_obj)));
        %initialObject = ifftshift(kernel).*(obj_temp - initialObject);
        %initialObject = obj_temp - initialObject;
    end
    
    u = initialObject + ds*(u - u_K);    
    
end

reconstructionTime = toc;
reconstructionTime = round(10*reconstructionTime)./10;
fprintf('GENFIRE: Reconstruction completed in %.12g seconds.\n\n',reconstructionTime);

%obj.reconstruction = My_stripzero(fftshift(obj.reconstruction),[obj.Dim1 obj.Dim2 obj.Dim1]);
%obj.final_rec = fftshift(initialObject);
%obj.final_rec = obj.final_rec(obj.Vol_ind(1,1):obj.Vol_ind(1,2), obj.Vol_ind(2,1):obj.Vol_ind(2,2), obj.Vol_ind(3,1):obj.Vol_ind(3,2));
end

%{
function obj = get_R(obj)
for projNum = 1:size(kMeasured,3);
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
    
    measuredX(:,:,projNum) = rotkCoords(1,:);%rotated X
    measuredY(:,:,projNum) = rotkCoords(2,:);%rotated Y
    measuredZ(:,:,projNum) = rotkCoords(3,:);%rotated Z
end
end
%}