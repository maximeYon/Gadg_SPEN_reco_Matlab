%% Reconstruction of SPEN data in .h5 generated with the EPI_SPEN export
clearvars; close all; clc;
addpath('/usr/local/share/ismrmrd/matlab')

%% Load the SPEN data
% filename='/home/mygadg/Data/RAW_EPI/meas_MID00106_FID14252_SPEN_diff_post_1seg.h5';
% filename='/home/mygadg/Data/RAW_EPI/meas_MID00108_FID14254_SPEN_diff_post_3seg.h5';
filename='/home/mygadg/Data/RAW_EPI/meas_MID00110_FID14256_SPEN_diff_post_4seg.h5';

dset = ismrmrd.Dataset(filename, 'dataset');
% D_noise = dset_noise.readAcquisition();
D = dset.readAcquisition();

%% Remove noise scans
isNoise = D.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
firstScan = find(isNoise==0,1,'first');
if firstScan > 1
    noise = D.select(1:firstScan-1);
else
    noise = [];
end

meas  = D.select(firstScan:D.getNumber);

clearvars D firstScan isNoise;

%% Get general parameters
% hdr_noise = ismrmrd.xml.deserialize(dset_noise.readxml);
hdr = ismrmrd.xml.deserialize(dset.readxml);
clearvars dset;
parameters = struct;

% Matrix size
parameters.Nro = hdr.encoding.encodedSpace.matrixSize.x;
parameters.Npe = hdr.encoding.encodedSpace.matrixSize.y;
parameters.enc_Z = hdr.encoding.encodedSpace.matrixSize.z;

% Field of View
parameters.dFOVReadout = hdr.encoding.encodedSpace.fieldOfView_mm.x;
parameters.dFOVPhase = hdr.encoding.encodedSpace.fieldOfView_mm.y;
parameters.dThickness = hdr.encoding.encodedSpace.fieldOfView_mm.z;

[ parameters.nSlices, parameters.number_of_channels  , parameters.number_of_repetitions, parameters.number_of_contrasts, parameters.number_of_phase, parameters.number_of_average , parameters.number_of_segments, parameters.number_of_sets] = get_number_of( hdr );

%% get flags
Flags = struct;
Flags.ACQ_IS_DUMMYSCAN_DATA = find(meas.head.flagIsSet('ACQ_IS_DUMMYSCAN_DATA')) ;
Flags.ACQ_FIRST_IN_SEGMENT = find(meas.head.flagIsSet('ACQ_FIRST_IN_SEGMENT')) ;
Flags.ACQ_FIRST_IN_SET = find(meas.head.flagIsSet('ACQ_FIRST_IN_SET')) ;
Flags.ACQ_IS_REVERSE = find(meas.head.flagIsSet('ACQ_IS_REVERSE')) ;
Flags.ACQ_IS_NOT_REVERSE = find(~meas.head.flagIsSet('ACQ_IS_REVERSE')) ;
Flags.ACQ_IS_PARALLEL_CALIBRATION = find(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')) ;
Flags.ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING = find(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING')) ;
Flags.ACQ_IS_NOT_PARALLEL_CALIBRATION = find(~meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')) ;
Flags.ACQ_IS_NAVIGATION_DATA = find(meas.head.flagIsSet('ACQ_IS_NAVIGATION_DATA')) ;
Flags.ACQ_IS_PHASECORR_DATA = find(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ;
Flags.all_ACQS = 1:size(meas.data,2);
Flags.Not_ACQS = sort([Flags.ACQ_IS_PHASECORR_DATA Flags.ACQ_IS_NAVIGATION_DATA]);
Flags.all_ACQS(Flags.Not_ACQS)=[];

% figure(1)
% subplot(311); plot(meas.head.idx.kspace_encode_step_1);
% % subplot(312);plot(meas.head.idx.kspace_encode_step_2);
% subplot(312);plot(meas.head.idx.set);
% subplot(313);plot(meas.head.idx.segment);

%% Get specific SPEN parameters
parameters.EchoSpacing = hdr.sequenceParameters.echo_spacing*1000;
parameters.LPE = parameters.dFOVReadout/10; % FOV along PE [cm]
parameters.nReps = parameters.number_of_repetitions;
Shifts = meas.head.position(:,Flags.all_ACQS(1));
parameters.PEShift = Shifts(2,1);

% get user parameters
for ind = 1:size(hdr.userParameters.userParameterLong,2)
    parameters.(hdr.userParameters.userParameterLong(1,ind).name) = hdr.userParameters.userParameterLong(1,ind).value;
end
clearvars ind Shifts;

% prendre slice 1
% acqs_test = find(meas.head.idx.slice==0 ) ;
% size(acqs_test)
% meas.head.position(3,acqs_test)  % aSlice.position.dTra
% 
% meas.head.read_dir 
% meas.head.phase_dir     % aSlice.position  Pas sur ?
% meas.head.slice_dir

%% test
% set = double(meas.head.idx.set(Flags.all_ACQS)+1);
% repetition = double(meas.head.idx.repetition(Flags.all_ACQS)+1);
% e2 = double(meas.head.idx.kspace_encode_step_1(Flags.all_ACQS)+1);
% n_enc = max(e2);
% nseg = max(set)/max(repetition);
% 
% figure()
% hold on
% plot(set,'r')
% plot(repetition ,'b')
% 
% 
% figure()
% plot(e2+((set-1)*n_enc-(repetition-1)*(n_enc*nseg)))



%% Filter meas to get complex data of measurement only
% clear meas and create Cmplxdata :
% dim1,kspace_encode_step_1,kspace_encode_step_2,slice,set(segment+diff),coils
CmplxData = zeros(size(meas.data{Flags.all_ACQS(1)},1),max(meas.head.idx.kspace_encode_step_1(Flags.all_ACQS))+1,max(meas.head.idx.kspace_encode_step_2(Flags.all_ACQS))+1,max(meas.head.idx.slice(Flags.all_ACQS))+1,max(meas.head.idx.set(Flags.all_ACQS))+1,parameters.number_of_channels);
for ind=Flags.all_ACQS
    CmplxData(:,meas.head.idx.kspace_encode_step_1(ind)+1,meas.head.idx.kspace_encode_step_2(ind)+1,meas.head.idx.slice(ind)+1,meas.head.idx.set(ind)+1,:) = meas.data{ind};
end

% retieve NShots and re arrange the data
parameters.nShots = double((max(meas.head.idx.set(Flags.all_ACQS))+1) / (max(meas.head.idx.repetition(Flags.all_ACQS))+1));
parameters.nWeight = double(max(meas.head.idx.repetition(Flags.all_ACQS))+1);

CmplxData = reshape(CmplxData,size(CmplxData,1),size(CmplxData,2),size(CmplxData,3),size(CmplxData,4),parameters.nShots,size(CmplxData,5)/parameters.nShots,size(CmplxData,6));
CmplxData = permute(CmplxData, [1 2 5 3 4 6 7]);
shot_coord = 1:parameters.nShots:parameters.nShots*size(CmplxData,2);
for shot = 1:parameters.nShots
CmplxDataSeg(:,shot_coord+shot-1,:,:,:,:) = CmplxData(:,:,shot,:,:,:,:);
end
CmplxDataSeg = reshape(CmplxDataSeg,size(CmplxData,1),size(CmplxDataSeg,2),size(CmplxData,4),size(CmplxData,5),size(CmplxData,6),size(CmplxData,7));
CmplxData = CmplxDataSeg;
clearvars CmplxDataSeg shot shot_coord ind;

CmplxData = permute(CmplxData,[1 2 6 4 3 5]);
CmplxData = reshape(CmplxData,size(CmplxData,1),size(CmplxData,2),size(CmplxData,3),size(CmplxData,4),size(CmplxData,5),1,size(CmplxData,6));

%% Remove oversampling
CmplxData = CmplxData(size(CmplxData,1)/4:size(CmplxData,1)-size(CmplxData,1)/4-1,:,:,:,:,:,:);
%% SPEN reconstruction
addpath(genpath('/home/mygadg/Documents/MATLAB/SPEN_Gadgetron_Reco/SPEN_Siemens_Recon/'))
% The Readout reggriding need to be performed by the gadgetron.

% CmplxData size(120 120 12 3 1 1 13)
% CmplxData = zeros(120,120,12,3,1,1,13);

%% user parameter
TVW_L2.Final=0.1;  % Weight for L2 regularization for final reconstruction
TVW_L2.EO=0.6;  % Weight for L2 regularization for reconstruction of even and odd images used to determine EO phase correction
TVW_L2.ImgPerShot=0.6;  % Weight for L2 regularization for reconstruction of per-shot images used to perform motion correction
TresholdCC=0.15; % Threshold to determine to how many coils to compress the data
WeightingFactor=100; % Weighting factor used for super-resolution (SR) matrix calculation

DoMotionCorrection=true; % Set true for in-vivo data, false for phantom data, relevant only for data acquired with #_of_segments > 1
[SPEN_Image]=SPEN_Reconstruction_Gadg(DoMotionCorrection,TVW_L2,TresholdCC,WeightingFactor,parameters,CmplxData);

% Display
figure()
for dim5 = 1: size(SPEN_Image,5)
    for dim3 = 1: size(SPEN_Image,3)
       subplot(1,size(SPEN_Image,3),dim3)
       imagesc(abs(squeeze(SPEN_Image(:,:,dim3,1,dim5))))
       colormap gray
    end
    pause(0.5)
end

