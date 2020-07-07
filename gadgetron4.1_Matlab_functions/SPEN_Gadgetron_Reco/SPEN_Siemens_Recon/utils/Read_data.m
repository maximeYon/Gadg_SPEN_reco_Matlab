% Script for reading data and getting some variables necessary for SR
% matrix calculation

PathData=[ComputerBasePath '/' DataFileName ];
disp(['Working on ' PathData]);
twix_obj=mapVBVD(PathData,'removeOS','ignoreSeg','rampSampRegrid');
image_data_unsorted = twix_obj{end}.image.unsorted();
sizeUnsorted=size(image_data_unsorted);


AccelFactor=twix_obj{end}.hdr.Phoenix.sPat.lAccelFactPE;
nSlices=twix_obj{end}.image.NSli;
nPE_acquired=length(1:AccelFactor:twix_obj{end}.image.NLin);


% get count of b0 and bweighted images
numberB0images=0;
if isfield(twix_obj{end}.hdr.Phoenix.sDiffusion,'alBValue')
    for i=1:size(twix_obj{end}.hdr.Phoenix.sDiffusion.alBValue,2) % 1:how many bvalues do we use
        if isempty(twix_obj{end}.hdr.Phoenix.sDiffusion.alBValue{i});
            nWeight{i}=1;
        else
            nWeight{i}=twix_obj{end}.hdr.Phoenix.sDiffusion.lDiffDirections;
        end
    end
    nWeight=sum(cell2mat(nWeight).*cell2mat(twix_obj{end}.hdr.Phoenix.sDiffusion.alAverages));
else
    numberB0images=1;
    numberBweightedImages=0;
    numberOfDiffDirections=0;
    nWeight=numberB0images;
end

if isfield(twix_obj{end}.hdr.Phoenix.sDiffusion,'alBValue')
    bvalues=twix_obj{end}.hdr.Phoenix.sDiffusion.alBValue;
end

if ~exist('bvalues')
    bvalues=[];
end

ProtShotsIdx=9;
nShots=twix_obj{end}.hdr.Phoenix.sWipMemBlock.alFree{ProtShotsIdx};
if isempty(nShots)
    nShots=twix_obj{end}.hdr.Phoenix.sWipMemBlock.alFree{7};
end

image_data_sorted_size= [sizeUnsorted(1) sizeUnsorted(2) nPE_acquired nSlices nShots nWeight];
nReps=sizeUnsorted(3)/prod(image_data_sorted_size(3:end));
image_data_sorted_size(end+1)=nReps;
image_data_sorted=reshape(image_data_unsorted, image_data_sorted_size); % [RO x nCoils x nPEacquired x nSlice x nShots x bWeight x nReps]

clear image_data_unsorted

BeforeCC=permute(image_data_sorted,[3 1 2 4 7 6 5]); % [nPEacquired x RO x nCoils x nSlice x nReps x bWeight x nShots]

clear image_data_sorted

Combined=CombineDims(BeforeCC(:,:,:,:,:,:,:),[1 7]);

clear BeforeCC
sizeCombined=size(Combined);
CmplxData=zeros([twix_obj{end}.image.NLin*nShots sizeCombined(2:end)]);
AcquiredLines=size(Combined,1):-AccelFactor:1;
CmplxData(AcquiredLines(end):AccelFactor:end,:,:,:,:,:,:)=Combined;

clear Combined

slicePos=twix_obj{end}.image.slicePos;
[~,sliceorder]=sort(unique(slicePos(3,:),'stable')); % check this: might change for ther orientations than transversal
if (length(sliceorder)==1 && nSlices~=1)
    [~,sliceorder]=sort(unique(slicePos(2,:),'stable')); % check this: might change for ther orientations than transversal
end

CmplxData=CmplxData(:,:,:,sliceorder,:,:,:,:);
CmplxData=permute(CmplxData,[1 2 3 4 5 7 6]);

% get parameters to determine PE shift
asSlice=twix_obj{end}.hdr.MeasYaps.sSliceArray.asSlice{1};

asSlice.sPosition
asSlice.sNormal

RotMat = transpose(Quat2RotMat(slicePos(4:end,1))) ;
asSlice=twix_obj{end}.hdr.MeasYaps.sSliceArray.asSlice{1};
RotMat = transpose(Quat2RotMat(slicePos(4:end,1))) ;
ImgCenterShifted=slicePos(1:3,1);
RotatedLocs=RotMat.'*ImgCenterShifted;
PEShift=RotatedLocs(1,1);
PEShiftSign=-1;
PEShift=PEShiftSign*PEShift;
rvalue=twix_obj{end}.hdr.Phoenix.sWipMemBlock.alFree{1}; % Rvalue of CHIRP
EchoSpacing=twix_obj{end}.hdr.Config.EchoSpacing_eff_us;% useconds
LPE=twix_obj{end}.hdr.Phoenix.sSliceArray.asSlice{1}.dPhaseFOV/10 ; % FOV along PE [cm]
Npe = size(CmplxData,1); % points along PE
Nro = size(CmplxData,2); % points along RO
nSlices=size(CmplxData,4); % number of slices