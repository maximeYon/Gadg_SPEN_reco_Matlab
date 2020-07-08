function [L2image]=FinalReconstruction(CmplxData, gccmtx_aligned_finalx, nccx, EvenOddMapStoredBySlices, FinalA, FinalAFull,...
    FinalB0Recon, nShots, TVW_L2, DoMotionCorrection,pFT)

% Function for reconstructing images either with L2 regularization.
% Depending if data is acquired with even or odd number of shots, splitting
% in even and odd subsets is slightly different, otherwise the recon is the same.
%
% The pipeline works like this:
% 1) Compress the data to one channel (to determine phase per shot for
% motion correction) or nccx channels (for final reconstruction)
% 2) Split data in even and odd datasets and apply even-odd correction
% along RO determined previously and stored in even-odd maps. Do this for
% data compressed to one and nccx channels.
% 3) Combine corrected even and odd data back together
% 4) Determine phase introduced by motion for each shot separately, first,
% by reconstructing image of final resolution for each shot separately
% using SUSPENSE principle.
% 5) Determine sensitivity maps from B0 data we reconstructed previously
% while determining even-odd correction.
% 6) Apply motion correction to sensitivity maps
% 7) Do iterative image reconstruction with L2  regularization


CCCmplxDataSx=ApplySCCBySlices(CmplxData(:,:,:,:,:,:,:),gccmtx_aligned_finalx,nccx); % Compressing to nccx channels
CCCmplxDataS1=ApplySCCBySlices(CmplxData(:,:,:,:,:,:,:),gccmtx_aligned_finalx,1); % compressing data to one channel

Nro=size(CmplxData,2); % number of lines in read-out direction
Npe=size(CmplxData,1); % number of lines in Phase-encoding direction
nSlices=size(CmplxData,4); % # of slices
nWeight=size(CmplxData,7); % # of bvalues
nReps=size(CmplxData,5); % # of repetitions

if pFT<1
%    CCCmplxDataSx(:,1(Nro*pFT+1):end,:,:,:,:,:)=0;
    CCCmplxDataSx(:,1:ceil((1-pFT)*Nro),:,:,:,:,:)=0;
    for SliI=1:nSlices
        for bval=1:nWeight
            for i=1:nccx
        [~, CCCmplxDataSx(:,:,i,SliI,:,bval)] = pocs(CCCmplxDataSx(:,:,i,SliI,:,bval), 20, 0);
            end
        end
    end
%    CCCmplxDataS1(:,(Nro*pFT+1):end,:,:,:,:,:)=0;
    CCCmplxDataS1(:,1:ceil((1-pFT)*Nro),:,:,:,:,:)=0;
    for SliI=1:nSlices
        for bval=1:nWeight
    [~, CCCmplxDataS1(:,:,:,SliI,:,bval)] = pocs(CCCmplxDataS1(:,:,:,SliI,:,bval), 20, 0);
        end
    end
end

L2image=zeros(Npe,size(CCCmplxDataS1,2),nSlices,nWeight,nReps);

if mod(nShots,2)==1
    if nSlices>3
        parfor_arg =Inf;
    else
        parfor_arg = 0;
    end
    parfor (SliI=1:nSlices, parfor_arg)
        
        for Ww=1:nWeight
            for RIndx=1:nReps
                CCCurCmplxDataSx=CCCmplxDataSx(:,:,:,SliI,RIndx,Ww); % select just one slice, all channels
                CCCurCmplxDataS1=CCCmplxDataS1(:,:,:,SliI,RIndx,Ww); % select one slice, one channel
                
                %% Apply EO
                % First we apply even-odd correction to data compressed to
                % nccx channels
                CCCurDataSE=CCCurCmplxDataSx(1:2:end,:,:,:,:,:,:,:); % even data, one slice
                CCCurDataSO=CCCurCmplxDataSx(2:2:end,:,:,:,:,:,:,:); % odd data, one slice
                
                CCCurDataSCOC=RepDotMult( fft1cg(CCCurDataSO,2),EvenOddMapStoredBySlices(1,:,SliI,:,:,:,:,:,:,:,:)); % Apply even-odd correction
                
                % Combine even and odd data data with EO correction applied
                Cat=zeros(size(CCCurCmplxDataSx));
                Cat(1:2:end,:,:,:,:,:,:)=CCCurDataSE;
                Cat(2:2:end,:,:,:,:,:,:)=ifft1cg(CCCurDataSCOC,2);
                
                CmplxDataEOx=Cat;
                Cat=[];
                
                % Apply even-odd correction to data compressed to one
                % channels.
                
                CCCurDataSE=CCCurCmplxDataS1(1:2:end,:,:,:,:,:,:,:); % get even data
                CCCurDataSO=CCCurCmplxDataS1(2:2:end,:,:,:,:,:,:,:); % get odd data
                
                CCCurDataSCOC=RepDotMult( fft1cg(CCCurDataSO,2),EvenOddMapStoredBySlices(1,:,SliI,:,:,:,:,:,:,:,:)); % Apply even odd
                Cat=zeros(size(CCCurCmplxDataS1));
                
                % Combine corrected even and odd data
                Cat(1:2:end,:,:,:,:,:,:,:)=CCCurDataSE;
                Cat(2:2:end,:,:,:,:,:,:,:)=ifft1cg(CCCurDataSCOC,2);
                
                CmplxDataEO1=Cat;
                Cat=[];
                %% Phase per shot
                
                % Determine phase introduced by motion for each shot
                
                SensShotCorrection=ones([size(FinalA,2) gsize(CCCurDataSE,2) nShots]);
                
                Cat=zeros([size(CmplxDataEO1,1)/nShots size(CmplxDataEO1,2) nShots]);
                % Splitting the data in separate shots
                for i=1:nShots
                    Cat(:,:,i)=CmplxDataEO1(i:nShots:end,:);
                end
                
                CurData=fft1cg(Cat,2);
                ImgPerShot=[];
                
                % Reconstructing high resolution image for each shot
                % separately
                for ShotI=1:nShots;
                    ImgPerShot(:,:,ShotI)=RunfnlCgIterationsf(1000+1,FinalA(:,:,ShotI),11,SensShotCorrection(:,:,1),TVW_L2.ImgPerShot,CurData(:,:,ShotI),2,struct('ShowFig',false));
                end
                
                
                PhaseDiff=exp(+1i*angle(ImgPerShot));
                
                %% Getting Sens Maps from data
                
                if nccx==1
                    SensE=ones(size(FinalB0Recon(:,:,:,SliI)));
                else
                    SensE= RunESPIRiTForSensMaps( FinalB0Recon(:,:,:,SliI),25,[6 6]);
                end
                
                %% Final Reconstruction
                
                CurData=fft1cg(PartitionDimInterleaved(CmplxDataEOx,1,nShots),2);
                
                % We apply previously determined phase introduced by motion to sensitivity maps
                
                CurSens=repmat(SensE,[1 1 1 1 nShots]); % sensitivity maps are assumed to be the same for each shot and repetition
                RPhaseDiff=permute(PhaseDiff,[1 2 5 4 3]);
                RPhaseDiff=repmat(RPhaseDiff,[1 1 nccx 1 1]); % Motion is the same for each channel
                % To do or not to do motion correction (for phantom will probably introduce artifacts)
                if DoMotionCorrection
                    CurSensC=CurSens.*RPhaseDiff; % Correcting sensitivity maps for motion
                else
                    CurSensC=CurSens;
                end
                % Preparing Super-resolution matrix needed as input for final reconstruction
                CurSR=repmat(FinalAFull,[1 1 nccx 1]);
                CurSR=PartitionDimInterleaved(CurSR,1,nShots);
                
                % Some data reordering for input for iterative
                % reconstruction
                CurD=CombineDims(CombineDims(CurData,[4 5]),[3 4]);
                CurS=CombineDims(CombineDims(CurSensC,[4 5]),[3 4]);
                CurA=CombineDims(CombineDims(CurSR,[4 5]),[3 4]);
                
                % Reconstructing your image.
                L2image(:,:,SliI,Ww,RIndx)=RunfnlCgIterationsf(1000+1,CurA,2,CurS,TVW_L2.Final,CurD,2,struct('ShowFig',false));%L2 including Normalization (MATLAB)
                
            end
            
        end
        
    end
else
    if nSlices>3
        parfor_arg =Inf;
    else
        parfor_arg = 0;
    end
    parfor (SliI=1:nSlices, parfor_arg)
        
        for Ww=1:nWeight
            for RIndx=1:nReps
                CCCurCmplxDataSx=CCCmplxDataSx(:,:,:,SliI,RIndx,:,Ww); % select just one slice, all channels
                CCCurCmplxDataS1=CCCmplxDataS1(:,:,:,SliI,RIndx,:,Ww); % select one slice, one channel
                
                %% Apply EO
                % First we apply even-odd correction to data compressed to
                % nccx channels
                CCCurDataSSplit=PartitionDimInterleaved(CCCurCmplxDataSx,1,nShots); % split data between shots
                CCCurDataSSplitE=CCCurDataSSplit(1:2:end,:,:,:,:,:,:,:); % even data, one slice
                CCCurDataSSplitO=CCCurDataSSplit(2:2:end,:,:,:,:,:,:,:); % odd data, one slice
                
                CCCurDataSCOC=RepDotMult( fft1cg(CCCurDataSSplitO,2),EvenOddMapStoredBySlices(1,:,SliI,:,:,:,:,:,:,:,:)); % Apply even-odd correction
                
                % Combine even and odd data data with EO correction applied
                Cat=CombineDataDim(CCCurDataSSplitE,ifft1cg(CCCurDataSCOC,2));
                nDims=length(size(squeeze(Cat)));
                CmplxDataEOx=CombineDims(squeeze(Cat),[1 nDims]);
                
                
                % Apply even-odd correction to data compressed to one
                % channels.
                CCCurDataSSplit=PartitionDimInterleaved(CCCurCmplxDataS1,1,nShots); % split between shots
                CCCurDataSSplitE=CCCurDataSSplit(1:2:end,:,:,:,:,:,:,:); % get even data
                CCCurDataSSplitO=CCCurDataSSplit(2:2:end,:,:,:,:,:,:,:); % get odd data
                
                CCCurDataSCOC=RepDotMult( fft1cg(CCCurDataSSplitO,2),EvenOddMapStoredBySlices(1,:,SliI,:,:,:,:,:,:,:,:)); % Apply even odd
                Cat=CombineDataDim(CCCurDataSSplitE,ifft1cg(CCCurDataSCOC,2));
                nDims=length(size(squeeze(Cat)));
                CmplxDataEO1=CombineDims(squeeze(Cat),[1 nDims]);
                
                %% Phase per shot
                
                % Determine phase introduced by motion for each shot
                
                SensShotCorrection=ones([size(FinalA,2) gsize(CCCurDataSSplitE,[2 3])]);
                CurData=permute(squeeze(fft1cg(Cat,2)),[1 2 4 3]);
                ImgPerShot=[];
                
                
                for ShotI=1:nShots;
                    ImgPerShot(:,:,ShotI)=RunfnlCgIterationsf(1000+1,FinalA(:,:,ShotI),11,SensShotCorrection(:,:,1),TVW_L2.ImgPerShot,CurData(:,:,ShotI),2,struct('ShowFig',false));
                end
                
                PhaseDiff=exp(+1i*angle(ImgPerShot));
                
                %% Getting Sens Maps from data
                
                if nccx==1
                    SensE=ones(size(FinalB0Recon(:,:,:,SliI)));
                else
                    SensE= RunESPIRiTForSensMaps( FinalB0Recon(:,:,:,SliI),25,[6 6]);
                end
                
                %% Final Reconstruction
                
                CurData=fft1cg(PartitionDimInterleaved(CmplxDataEOx,1,nShots),2);
                
                % We apply previously determined phase introduced by motion to sensitivity maps
                
                CurSens=repmat(SensE,[1 1 1 1 nShots]);
                RPhaseDiff=permute(PhaseDiff,[1 2 5 4 3]);
                RPhaseDiff=repmat(RPhaseDiff,[1 1 nccx 1 1]);
                
                % To do or not to do motion correction (for phantom will probably introduce artifacts)
                if DoMotionCorrection
                    CurSensC=CurSens.*RPhaseDiff; % Correcting sensitivity maps for motion
                else
                    CurSensC=CurSens;
                end
                
                % Preparing Super-resolution matrix needed as input for final reconstruction
                
                CurSR=repmat(FinalAFull,[1 1 nccx 1]);
                CurSR=PartitionDimInterleaved(CurSR,1,nShots);
                
                % Some data reordering for input for iterative
                % reconstruction
                
                CurD=CombineDims(CombineDims(CurData,[4 5]),[3 4]);
                CurS=CombineDims(CombineDims(CurSensC,[4 5]),[3 4]);
                CurA=CombineDims(CombineDims(CurSR,[4 5]),[3 4]);
                
                
                % Reconstructing your image.
                L2image(:,:,SliI,Ww,RIndx)=RunfnlCgIterationsf(1000+1,CurA(:,:,:),2,CurS(:,:,:),TVW_L2.Final,CurD(:,:,:),2,struct('ShowFig',false));%L2 including Normalization (MATLAB)
                
            end
        end
    end
end