function [Parameters] = reconstruct_SPEN_thermo(input,connection,Parameters)
%image is a structure navigating throug the functions
image = input(); % Call input function to produce the next slice.
acq_header = connection.header;
%% Define a repetition counter
counter = max(unique(image.bits.buffer.headers.set));

% Get the complex Data, no need to store them
CmplxData = image.bits.buffer.data;

tic
%% Retrieve the parameters and compute the super resolution matrix
%% used for deconvolution of all sets
if counter==0
    addpath('/usr/local/share/ismrmrd/matlab')
    addpath(genpath('/home/imagerie/Dev/Gadg_SPEN_reco_Matlab/gadgetron4.1_Matlab_functions/SPEN_Gadgetron_Reco'));
    
    %% Retrieve common parameters and store them in acq_header. structure
    
    acq_header.SPEN_parameters.matrixSize(1,1) = acq_header.encoding.encodedSpace.matrixSize.x;
    acq_header.SPEN_parameters.matrixSize(1,3) = acq_header.encoding.encodedSpace.matrixSize.z;
    acq_header.SPEN_parameters.FOV(1,1) = acq_header.encoding.encodedSpace.fieldOfView_mm.x;
    acq_header.SPEN_parameters.FOV(1,2) = acq_header.encoding.encodedSpace.fieldOfView_mm.y;
    acq_header.SPEN_parameters.FOV(1,3) = acq_header.encoding.encodedSpace.fieldOfView_mm.z;
    
    
    acq_header.SPEN_parameters.Enc1 = acq_header.encoding.encodingLimits.kspace_encoding_step_1.maximum+1;
    acq_header.SPEN_parameters.Enc2 = acq_header.encoding.encodingLimits.kspace_encoding_step_2.maximum+1;
    acq_header.SPEN_parameters.nSlices = acq_header.encoding.encodingLimits.slice.maximum+1;
    acq_header.SPEN_parameters.set = acq_header.encoding.encodingLimits.set.maximum+1;
    acq_header.SPEN_parameters.EchoSpacing = acq_header.sequenceParameters.echo_spacing*1000;
    
    % get user parameters
    for ind = 1:size(acq_header.userParameters.userParameterLong,2)
        acq_header.SPEN_parameters.(acq_header.userParameters.userParameterLong(1,ind).name) = acq_header.userParameters.userParameterLong(1,ind).value;
    end
    clearvars ind
    
    acq_header.SPEN_parameters.matrixSize(1,2) = acq_header.SPEN_parameters.Nky;
    %     acq_header.SPEN_parameters.nShots = acq_header.SPEN_parameters.Nseg;
    acq_header.SPEN_parameters.nShots = acq_header.SPEN_parameters.Nky/acq_header.SPEN_parameters.Enc1;
    acq_header.SPEN_parameters.repetition = acq_header.SPEN_parameters.set/acq_header.SPEN_parameters.nShots;
    if isfield(acq_header.SPEN_parameters,'RO_pFT')==1
        acq_header.SPEN_parameters.RO_pFT = acq_header.SPEN_parameters.RO_pFT./100;
    else
        acq_header.SPEN_parameters.RO_pFT =1;
    end
    
    acq_header.SPEN_parameters.bvalues = [acq_header.SPEN_parameters.bvalue_1  acq_header.SPEN_parameters.bvalue_2 acq_header.SPEN_parameters.bvalue_3  acq_header.SPEN_parameters.bvalue_4 acq_header.SPEN_parameters.bvalue_5  acq_header.SPEN_parameters.bvalue_6];
    acq_header.SPEN_parameters.bvalues = acq_header.SPEN_parameters.bvalues(1,1:acq_header.SPEN_parameters.DiffWeightings);
    
    %% Adapt specific SPEN parameters
    acq_header.SPEN_parameters.LPE = acq_header.SPEN_parameters.FOV(1,2)/10;
    acq_header.SPEN_parameters.nReps = 1; %% Need to be changed !!!!!!!!!!!!!!!!!
    acq_header.SPEN_parameters.Npe = acq_header.SPEN_parameters.matrixSize(1,2);
    acq_header.SPEN_parameters.nWeight = acq_header.SPEN_parameters.repetition;
    
    %% Define user parameter
    acq_header.SPEN_parameters.TVW_L2.Final=0.15;  % Weight for L2 regularization for final reconstruction
    acq_header.SPEN_parameters.TVW_L2.EO=0.6;  % Weight for L2 regularization for reconstruction of even and odd images used to determine EO phase correction
    acq_header.SPEN_parameters.TVW_L2.ImgPerShot=0.6;  % Weight for L2 regularization for reconstruction of per-shot images used to perform motion correction
    acq_header.SPEN_parameters.TresholdCC=0.15; % Threshold to determine to how many coils to compress the data
    acq_header.SPEN_parameters.WeightingFactor=100; % Weighting factor used for super-resolution (SR) matrix calculation
    acq_header.SPEN_parameters.DoMotionCorrection=true; % Set true for in-vivo data, false for phantom data, relevant only for data acquired with #_of_segments > 1
    
    %% Get parameters to compute PEshift
    % PE shift is independant to the slice position
    idx_data=find(image.bits.buffer.headers.kspace_encode_step_1 ~=0);
    idx_data = idx_data(1,1);
    
    % Get the position in read phase slice
    [~,idx_2,idx_3,idx_4,idx_5,idx_6] = ind2sub(size(image.bits.buffer.headers.kspace_encode_step_1),idx_data);
    PosXYZ = image.bits.buffer.headers.position(:,idx_2,idx_3,idx_4,idx_5,idx_6);
    
    % Get the orientation
    RotMatread = image.bits.buffer.headers.read_dir(:,idx_2,idx_3,idx_4,idx_5,idx_6);
    RotMatphase = image.bits.buffer.headers.phase_dir(:,idx_2,idx_3,idx_4,idx_5,idx_6);
    RotMatslice = image.bits.buffer.headers.slice_dir(:,idx_2,idx_3,idx_4,idx_5,idx_6);
    RotMat = [RotMatphase,RotMatread,RotMatslice]; % Is is order alway good ?
    clearvars idx_2 idx_3 idx_4 idx_5 idx_6 idx_data;
    
    % Compute PEShift
    RotatedLocs=RotMat.'*PosXYZ;
    PEShift=RotatedLocs(1,1);
    PEShiftSign=-1;
    acq_header.SPEN_parameters.PEShift=PEShiftSign*PEShift;
    acq_header.SPEN_parameters.RotMat = RotMat; % store RotMat used to compute slice position
    clearvars RotMatread RotMatphase RotMatslice RotMat RotatedLocs PEShift PEShiftSign PosXYZ;
    
    %% Assign the parameters used for the Super-Resolution matrix
    acq_header.ParamSR.LPE=acq_header.SPEN_parameters.LPE; % FOV along PE direction [cm]
    acq_header.ParamSR.rvalue=-acq_header.SPEN_parameters.rvalue; % CHIRP rvalue
    acq_header.ParamSR.EchoSpacing=acq_header.SPEN_parameters.EchoSpacing; % Echo spacing [us]
    acq_header.ParamSR.PEShift=acq_header.SPEN_parameters.PEShift; % Shift for data along PE direction to remove shift added by EPI read-out
    
    
    %% SuperResolution (SR) matrix calculation
    % SR matrix depends on the sequence in use. This
    % calculation is for 180deg encoding CHIRP and another 180deg
    % slice-non-selective inversion pulse. If you use different sequence, you might need
    % to change some parameters in the function (or write your own calculation
    % of SR matrix).
    acq_header.SPEN_parameters.ChirpAmpExtraFactor=1;
    
    % this is for the final reconstruction (size equal to final data size)
    [acq_header.FinalAFull]=SuperResolution_VariableSize_Sim_V3(acq_header.SPEN_parameters.Npe,acq_header.ParamSR,1,acq_header.SPEN_parameters.WeightingFactor,acq_header.SPEN_parameters.nShots,acq_header.SPEN_parameters.ChirpAmpExtraFactor);
    
    %this is for phase per shot used for motion correction
    acq_header.FinalA=PartitionDimInterleaved(acq_header.FinalAFull,1,acq_header.SPEN_parameters.nShots);
    
    Parameters.ParamSR = acq_header.ParamSR;
    Parameters.FinalAFull = acq_header.FinalAFull;
    Parameters.FinalA = acq_header.FinalA;
else
    acq_header.ParamSR = Parameters.ParamSR;
    acq_header.SPEN_parameters = Parameters.SPEN_parameters;
    acq_header.FinalAFull = Parameters.FinalAFull;
    acq_header.FinalA = Parameters.FinalA;
end

%% Perform the reconstruction for each set of image
counter = counter+1;
disp([' number of set = ' num2str(counter)])
%% Create CmplxData size(120 120 12 3 1 1 13) for SPENreco
% Crop the data in dim 2 : only the center contain the data due
% to the segmented acquistion and wrong tagging

CmplxData = CmplxData(:,size(CmplxData,2)/2-acq_header.SPEN_parameters.Enc1/2+1:size(CmplxData,2)/2+acq_header.SPEN_parameters.Enc1/2,:,:,:,:,:,:);
% Rearrange the segment
CmplxData = permute(CmplxData, [1 2 5 3 4 6 7]);
% %% Re-order the segments if sorting_dimension in AccTrig left blank
% acq_header.SPEN_parameters.SegOrder = unique(image.bits.buffer.headers.repetition)+1;
% CmplxData(:,:,1:max(acq_header.SPEN_parameters.SegOrder),:,:,:,:,:)= CmplxData(:,:,acq_header.SPEN_parameters.SegOrder,:,:,:,:,:);

shot_coord = 1:acq_header.SPEN_parameters.nShots:acq_header.SPEN_parameters.nShots*size(CmplxData,2);
CmplxDataSeg = zeros(size(CmplxData,1),size(CmplxData,2)*size(CmplxData,3),size(CmplxData,4),size(CmplxData,5),size(CmplxData,6),size(CmplxData,7));
for shot = 1:acq_header.SPEN_parameters.nShots
    CmplxDataSeg(:,shot_coord+shot-1,:,:,:,:) = CmplxData(:,:,shot,:,:,:,:);
end
CmplxDataSeg = reshape(CmplxDataSeg,size(CmplxData,1),size(CmplxDataSeg,2),size(CmplxData,5),size(CmplxData,7),1,1,size(CmplxData,8));
CmplxData = CmplxDataSeg;
clearvars CmplxDataSeg shot shot_coord;

%% permute dim 1 and  2
CmplxData = permute(CmplxData,[2 1 3 4 5 6 7]);

%% SPEN reconstruction
% PE Shift: Shift data along PE direction to remove phase added by EPI read-out
CmplxData = PEShiftingDataset(CmplxData,-acq_header.SPEN_parameters.PEShift,acq_header.SPEN_parameters.LPE,acq_header.SPEN_parameters.nShots);

%% pFT test and fill - mo
if acq_header.SPEN_parameters.RO_pFT<1.0
    CmplxData=pFT_Pad(CmplxData,acq_header.SPEN_parameters.RO_pFT);
end

%% Coil Compression
% Determine number of virtual channels to which compress data, but still
% not  to loose too much information (TresholdCC). If this condition,
% is not reached, we will compress to three coils.
if counter==1
    [~, acq_header.SPEN_parameters.gccmtx_aligned_finalx, Energyx]=CoilCompression(permute(CmplxData,[2 1 3 4 5 6 7]),1);
    acq_header.SPEN_parameters.nccx=find(Energyx(:,1)<acq_header.SPEN_parameters.TresholdCC,1);
    if isempty(acq_header.SPEN_parameters.nccx) && size(CmplxData,3)>2
        acq_header.SPEN_parameters.nccx=3;
    elseif isempty(acq_header.SPEN_parameters.nccx)
        acq_header.SPEN_parameters.nccx=1;
    end
    clearvars Energyx;
end

%% ONE SHOT Reconstruction
% Description of the workflow given inside the script. You might want to
% change L2 weight.
if acq_header.SPEN_parameters.nShots==1
    [SPEN_Image]=function_OneShotReconstruction(CmplxData,acq_header.SPEN_parameters.gccmtx_aligned_finalx,acq_header.SPEN_parameters.nccx,acq_header.FinalAFull,acq_header.SPEN_parameters,acq_header.SPEN_parameters.RO_pFT);
    % OneShotReconstruction
else
    %% Determine SENS map, as well as EO map
    % Determine EO maps for each slice from B0 data
    % The even-odd maps are assumed not to change for other b-weigted images and repetitions,
    % hence, they will be later applied for them also.
    % More explanations inside the function
    if counter==1
        acq_header.SPEN_parameters.RepsToDoEOAndB0=acq_header.SPEN_parameters.nReps;
        [acq_header.SPEN_parameters.EvenOddMapStoredBySlices, acq_header.SPEN_parameters.FinalB0Recon]=SensAndEOMapFromB0(CmplxData, acq_header.FinalAFull, acq_header.SPEN_parameters.gccmtx_aligned_finalx, acq_header.SPEN_parameters.nccx, acq_header.SPEN_parameters.nShots,acq_header.SPEN_parameters.RepsToDoEOAndB0,acq_header.SPEN_parameters.TVW_L2);
    end
    %% Final reconstruction
    [SPEN_Image]=FinalReconstruction(CmplxData, acq_header.SPEN_parameters.gccmtx_aligned_finalx, acq_header.SPEN_parameters.nccx, acq_header.SPEN_parameters.EvenOddMapStoredBySlices, acq_header.FinalA, acq_header.FinalAFull, acq_header.SPEN_parameters.FinalB0Recon, acq_header.SPEN_parameters.nShots, acq_header.SPEN_parameters.TVW_L2, acq_header.SPEN_parameters.DoMotionCorrection,acq_header.SPEN_parameters.RO_pFT);
end

SPEN_Image=permute(SPEN_Image, [2,1,3]);
SPEN_Image = flip(flip(SPEN_Image,2),1);

%% Siemens love square pixels, we need to regrid the image to fullfill this condition
size_x=acq_header.SPEN_parameters.FOV(1,1) / size(SPEN_Image,1);
size_y=acq_header.SPEN_parameters.FOV(1,2) / size(SPEN_Image,2);  
   
[X,Y] = meshgrid(1:size(SPEN_Image,1),1:size(SPEN_Image,2));
   
 if size_x>size_y
        fill_to=round((size_x/size_y)*size(SPEN_Image,1));
        [Xq,Yq] = meshgrid(1:(size(SPEN_Image,1)-1)/(fill_to-1):size(SPEN_Image,1),1:size(SPEN_Image,2));
        for i=1:size(SPEN_Image,3)
            SPEN_Image_temp(:,:,i) = interp2(X,Y,permute(SPEN_Image(:,:,i),[2 1 3]),Xq,Yq,'cubic');
        end
        SPEN_Image=permute(SPEN_Image_temp,[2 1 3]);
        clear SPEN_Image_temp
    else
        fill_to=round((size_y/size_x)*size(SPEN_Image,2));
        [Xq,Yq] = meshgrid(1:size(SPEN_Image,1),1:(size(SPEN_Image,2)-1)/(fill_to-1):size(SPEN_Image,2));
        for i=1:size(SPEN_Image,3)
            SPEN_Image_temp(:,:,i) = interp2(X,Y,permute(SPEN_Image(:,:,i),[2 1 3]),Xq,Yq,'cubic');
        end
        SPEN_Image=permute(SPEN_Image_temp,[2 1 3]);
        clear SPEN_Image_temp
    end

% Normalisation of the SPEN images (on the first one: b0 in diff)
if counter==1
    acq_header.SPEN_parameters.max_intensity=max(abs(SPEN_Image(:)));
end
SPEN_Image_no_norm = SPEN_Image;
SPEN_Image=uint16(4096*0.99*abs(SPEN_Image)/acq_header.SPEN_parameters.max_intensity);

SPEN_Image = single(SPEN_Image);
SPEN_Image = reshape(SPEN_Image,1,size(SPEN_Image,1),size(SPEN_Image,2),size(SPEN_Image,3));
SPEN_Image_no_norm = reshape(SPEN_Image_no_norm,1,size(SPEN_Image_no_norm,1),size(SPEN_Image_no_norm,2),size(SPEN_Image_no_norm,3));

time_reco = toc;
disp(['Duration of SPEN reco = ' num2str(time_reco) ' s'])

%% save the constant parameter to Parameter structure
Parameters.SPEN_parameters = acq_header.SPEN_parameters;

%% create suitable header for each slice
img_head = ismrmrd.ImageHeader;

%% Set the good header parameters for all slices
idx_data = find(image.bits.buffer.headers.kspace_encode_step_1 ~= 0);
idx_data = idx_data(1,1);
[~,idx_2,idx_3,idx_4,idx_5,idx_6] = ind2sub(size(image.bits.buffer.headers.kspace_encode_step_1),idx_data);

img_head.matrix_size(1) = size(SPEN_Image,2);
img_head.matrix_size(2) = size(SPEN_Image,3);
img_head.matrix_size(3) = size(SPEN_Image,4);

img_head.read_dir = image.bits.buffer.headers.read_dir(:,idx_2,idx_3,idx_4,idx_5,idx_6);
img_head.phase_dir = image.bits.buffer.headers.phase_dir(:,idx_2,idx_3,idx_4,idx_5,idx_6);
img_head.slice_dir = image.bits.buffer.headers.slice_dir(:,idx_2,idx_3,idx_4,idx_5,idx_6);
img_head.patient_table_position = image.bits.buffer.headers.patient_table_position(:,idx_2,idx_3,idx_4,idx_5,idx_6);
img_head.acquisition_time_stamp = image.bits.buffer.headers.acquisition_time_stamp(:,idx_2,idx_3,idx_4,idx_5,idx_6);
img_head.channels = 1;
img_head.data_type= 5;

%% retrieve slice order from positions
SliceShift = zeros(1,size(SPEN_Image,4));
for s = 1:size(SPEN_Image,4)
    idx_Encode = find(image.bits.buffer.headers.kspace_encode_step_1 ~= 0);
    idx_Slice = find(image.bits.buffer.headers.slice==s-1);
    idx_data = intersect(idx_Encode,idx_Slice);
    idx_data = idx_data(1,1);
    [~,idx_2,idx_3,idx_4,idx_5,idx_6] = ind2sub(size(image.bits.buffer.headers.kspace_encode_step_1),idx_data);
    Positions = image.bits.buffer.headers.position(:,idx_2,idx_3,idx_4,idx_5,idx_6);
    RotatedPositions=acq_header.SPEN_parameters.RotMat  .'*Positions;
    SliceShift(1,s)=RotatedPositions(3,1);
end
[~,Slice_order] = sort(SliceShift);

%% Set the good header parameters for each slice
slice_ind = 0;
for s = Slice_order
    idx_Encode = find(image.bits.buffer.headers.kspace_encode_step_1 ~= 0);
    idx_Slice = find(image.bits.buffer.headers.slice==s-1);
    idx_data = intersect(idx_Encode,idx_Slice);
    idx_data = idx_data(1,1);
    [~,idx_2,idx_3,idx_4,idx_5,idx_6] = ind2sub(size(image.bits.buffer.headers.kspace_encode_step_1),idx_data);
    img_head.position = image.bits.buffer.headers.position(:,idx_2,idx_3,idx_4,idx_5,idx_6);
    img_head.slice = image.bits.buffer.headers.slice(:,idx_2,idx_3,idx_4,idx_5,idx_6);
    % Unfortunatly this do not correct the slice position in case of
    % acquisition "entrelacé" (interleaved).
    
    %% Prepare image.data
    image_saved = gadgetron.types.Image.from_data(abs(SPEN_Image(:,:,:,s)),img_head);
    
    %% Prepare image.header
    image_saved.header.field_of_view(1,1) = acq_header.SPEN_parameters.FOV(1,1);
    image_saved.header.field_of_view(1,2) = acq_header.SPEN_parameters.FOV(1,2);
    image_saved.header.field_of_view(1,3) = acq_header.SPEN_parameters.FOV(1,3);
    
    image_saved.header.image_type = gadgetron.types.Image.MAGNITUDE;
    image_saved.header.image_series_index = 1000;
    slice_ind = slice_ind+1;
    image_saved.header.image_index = slice_ind+((counter-1)*Parameters.SPEN_parameters.nSlices)+1000;
    
    %% Send image
    connection.send(image_saved);
end

%% Send phase images
slice_ind = 0;
for s = Slice_order
    idx_Encode = find(image.bits.buffer.headers.kspace_encode_step_1 ~= 0);
    idx_Slice = find(image.bits.buffer.headers.slice==s-1);
    idx_data = intersect(idx_Encode,idx_Slice);
    idx_data = idx_data(1,1);
    [~,idx_2,idx_3,idx_4,idx_5,idx_6] = ind2sub(size(image.bits.buffer.headers.kspace_encode_step_1),idx_data);
    img_head.position = image.bits.buffer.headers.position(:,idx_2,idx_3,idx_4,idx_5,idx_6);
    img_head.slice = image.bits.buffer.headers.slice(:,idx_2,idx_3,idx_4,idx_5,idx_6);
    % Unfortunatly this do not correct the slice position in case of
    % acquisition "entrelacé" (interleaved).
    
    %% Prepare image.data
    image_saved = gadgetron.types.Image.from_data(single(angle(SPEN_Image_no_norm(:,:,:,s))),img_head);
    
    %% Prepare image.header
    image_saved.header.field_of_view(1,1) = acq_header.SPEN_parameters.FOV(1,1);
    image_saved.header.field_of_view(1,2) = acq_header.SPEN_parameters.FOV(1,2);
    image_saved.header.field_of_view(1,3) = acq_header.SPEN_parameters.FOV(1,3);
    
    image_saved.header.image_type = gadgetron.types.Image.PHASE;
    image_saved.header.image_series_index = 3000;
    slice_ind = slice_ind+1;
    image_saved.header.image_index = slice_ind+((counter-1)*Parameters.SPEN_parameters.nSlices)+3000;
    
    %% Send image
    connection.send(image_saved);
end

end

