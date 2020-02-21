function [EvenOddMapStoredBySlices, FinalB0Recon]=SensAndEOMapFromB0(CmplxDataSh, FinalAFull, gccmtx_aligned_finalx, nccx, nShots,nReps,TVW_L2)
% CmplxDataSh is the data in Fourier domaine
% FinalAFull should be nPE x nPE matrix (SUSPENSE EO)
% nReps is choosen in case one of the repetitions is corrupted, by default,
% last repetition is chosen
%
% The workflow of the function performing EO correction:
% 1) compress data to one channel;
%
% if nshots == odd
% 2) take each slice compressed to one channel and split  it in even and odd data sets;
% 3) split each even and odd data set in separate shots
% elseif nshots == even
%
% 2) split data set in separate shots
% 3) take each shot and split in even and odd data sets;
% end
%
% 4) split and repmat SR matrix, that each it's line correspond to the
% correct subset of spitted data;
% 5) Reconstruct two images from even and odd data separately with
% resolution of final image size;
% 6) determine phase difference between even and odd
% 7) do phase fitting to get linear phase difference along readout and
% construct phase maps;
% 8) Split data compressed to nccx coils in even and odd
% 9) Apply phase correction to this data
% 10) Reconstruct B0 image for each of the channels by simple
% multiplication with SR matrix and fft along RO direction: these are
% also our sensitivity maps

%%
nSlices=size(CmplxDataSh,4);
Npe=size(CmplxDataSh,1);
%% Coil Compressing
CCCurCmplxDataSx=ApplySCCBySlices(CmplxDataSh(:,:,:,:,1,1),gccmtx_aligned_finalx,nccx); % data compressed to nccx channels
CCCurCmplxDataS1=ApplySCCBySlices(CmplxDataSh(:,:,:,:,1,1),gccmtx_aligned_finalx,1); % data compressed to one channel


clear SensFull
clear FinalB0Recon
%%

if mod(nShots,2)==1
    if nSlices>2
        parfor_arg =0;
    else
        parfor_arg = Inf;
    end
    parfor (SliI=1:nSlices, parfor_arg)
        %% Even Odd Correction
        CurDataEO=CCCurCmplxDataS1(:,:,:,SliI,:,:,:);
        
        CCCurDataE=CurDataEO(1:2:end,:,:,:,:,:,:,:); % get even data
        CCCurDataO=CurDataEO(2:2:end,:,:,:,:,:,:,:); % get odd data
        
        CCCurDataSSplitE=zeros(size(CCCurDataE,1)/nShots,size(CCCurDataE,2),nShots);
        CCCurDataSSplitO=zeros(size(CCCurDataO,1)/nShots,size(CCCurDataO,2),nShots);
        
        for i=1:nShots
            CCCurDataSSplitE(:,:,i)=CCCurDataE(i:nShots:end,:); % split in shots even data
            CCCurDataSSplitO(:,:,i)=CCCurDataO(i:nShots:end,:); % split in shots odd data
        end
        
        FinalATryE=FinalAFull(1:2:end,:,:,:,:); % get SR matrix for even data
        FinalATryO=FinalAFull(2:2:end,:,:,:,:); % get SR matrix for odd data
        
        FinalATrySplitE=zeros(size(FinalATryE,1)/nShots,size(FinalATryE,2),nShots);
        FinalATrySplitO=zeros(size(FinalATryO,1)/nShots,size(FinalATryO,2),nShots);
        
        for i=1:nShots
            FinalATrySplitE(:,:,i)=FinalATryE(i:nShots:end,:); % split even SR matrix for each shot
            FinalATrySplitO(:,:,i)=FinalATryO(i:nShots:end,:); % split odd SR matrix for each shot
        end
        
        SensShotCorrection=ones([size(FinalATrySplitE,2) gsize(CCCurDataSSplitE,[2 3 4 5 6 7 8 9 10])]);
        
        % might use different TVW for regularization (0.1 is somewhat arbitrary choice, but seems to work fine)
        FinalE=RunfnlCgIterationsf(1000+1,FinalATrySplitE,11,squeeze(SensShotCorrection(:,:,:,:,:,:,:)),TVW_L2.EO,squeeze(fft1cg(CCCurDataSSplitE(:,:,:,:,:,:,:),2)),2,struct('ShowFig',false)); % Reconstruct Even image
        FinalO=RunfnlCgIterationsf(1000+1,FinalATrySplitO,11,squeeze(SensShotCorrection(:,:,:,:,:,:,:)),TVW_L2.EO,squeeze(fft1cg(CCCurDataSSplitO(:,:,:,:,:,:,:),2)),2,struct('ShowFig',false)); % Reconstruc Odd image
        
        %%
        CurW=abs(FinalE);
        SRPhaseDiff=angle(FinalE)-angle(FinalO); %phase difference between odd and even data
        EOWithW=CurW.*exp(1i*SRPhaseDiff); % some weight on the phase
        [ROEvenOddCorrection_bvalu, BestX, BestCost] =Robust1DLinearPhaseEstimate(EOWithW,2); % do a 1D linear fitting on the phase that will later be applied to all readout lines
        EvenOddMap=exp(+1i*angle(repmat(ROEvenOddCorrection_bvalu(1,:,:,:,:,:,:,:),[Npe 1 1 1 1 ])));
        EvenOddMapStoredBySlices(:,:,SliI)=EvenOddMap;
        
        CCCurDataE=[];
        CCCurDataO=[];
        CCCurDataSSplitE=[];
        CCCurDataSSplitO=[];
        
        % Apply EO to data compressed to nccx channels
        CCCurDataE=CCCurCmplxDataSx(1:2:end,:,:,SliI,:,:,:,:); % Even data
        CCCurDataO=CCCurCmplxDataSx(2:2:end,:,:,SliI,:,:,:,:); % Odd data
        
        CCCurDataSCOC=RepDotMult( fft1cg(CCCurDataO,2),EvenOddMap(1,:,:,:,:,:,:,:,:,:,:)); % applying even-Odd correction
        
        Cat=zeros([size(CCCurDataE,1)+size(CCCurDataO,1) gsize(CCCurDataE,2:4)]);
        Cat(1:2:end,:,:,:,:)=CCCurDataE;
        Cat(2:2:end,:,:,:,:)=ifft1cg(CCCurDataSCOC,2); % recombining the data
        
        DataFinal=fft1cg(Cat,2); % FT along RO
        
        SensFull=MultMatTensor(FinalAFull',DataFinal); % Multiplication with SR matrix to Even-Odd corrected data and
        % getting our image, that will serve as a sensitivity map at the same time for each channel
        
        
        FinalB0Recon(:,:,:,SliI)=SensFull;
    end
else % follow explanations above for odd number of shots: the same thing goes on there, just in the case of even shots,
    % the direction of readout is the same for n=nShots consecutive
    % lines, so splitting in even and odd data sets is slightly
    % different.
    
    if nSlices>2
        parfor_arg =0;
    else
        parfor_arg = Inf;
    end
    parfor (SliI=1:nSlices, parfor_arg)
        
        %% Even Odd Correction
        CurDataEO=CCCurCmplxDataS1(:,:,:,SliI,:,:);
        CCCurDataSSplit=PartitionDimInterleaved(CurDataEO,1,nShots);
        
        CCCurDataSSplitE=CCCurDataSSplit(1:2:end,:,:,:,:,:,:,:);
        CCCurDataSSplitO=CCCurDataSSplit(2:2:end,:,:,:,:,:,:,:);
        
        FinalATrySplit=PartitionDimInterleaved(FinalAFull,1,nShots);
        
        FinalATrySplitE=FinalATrySplit(1:2:end,:,:,:,:);
        FinalATrySplitO=FinalATrySplit(2:2:end,:,:,:,:);
        SensShotCorrection=ones([size(FinalATrySplitE,2) gsize(CCCurDataSSplitE,[2 3 4 5 6 7 8 9 10])]);
        
        % might use different TVW for regularization (0.1 is somewhat arbitrary choice, but seems to work fine)
        FinalE=RunfnlCgIterationsf(1000+1,FinalATrySplitE,11,squeeze(SensShotCorrection(:,:,:,:,:,:,:,:)),TVW_L2.EO,squeeze(fft1cg(CCCurDataSSplitE(:,:,:,:,:,:,:,:),2)),2,struct('ShowFig',false)); % Reconstruct Even image
        FinalO=RunfnlCgIterationsf(1000+1,FinalATrySplitO,11,squeeze(SensShotCorrection(:,:,:,:,:,:,:,:)),TVW_L2.EO,squeeze(fft1cg(CCCurDataSSplitO(:,:,:,:,:,:,:,:),2)),2,struct('ShowFig',false)); % Reconstruct Odd image
        %%
        CurW=abs(FinalE);
        SRPhaseDiff=angle(FinalE)-angle(FinalO); %pahse difference between odd and even
        EOWithW=CurW.*exp(1i*SRPhaseDiff); % some weight on the phase
        [ROEvenOddCorrection_bvalu, BestX, BestCost] =Robust1DLinearPhaseEstimate(EOWithW,2); % do a 1D linear fit on the phase that will later be applied to all readout lines
        EvenOddMap=exp(+1i*angle(repmat(ROEvenOddCorrection_bvalu(1,:,:,:,:,:,:,:),[Npe 1 1 1 1 ])));
        EvenOddMapStoredBySlices(:,:,SliI)=EvenOddMap;
        
        %% Apply EO
        CurDataChSh=PartitionDimInterleaved(CCCurCmplxDataSx(:,:,:,SliI,:,1,:),1,nShots);
        CCCurDataSSplitE=CurDataChSh(1:2:end,:,:,:,:,:,:,:);
        CCCurDataSSplitO=CurDataChSh(2:2:end,:,:,:,:,:,:,:);
        
        CCCurDataSCOC=RepDotMult( fft1cg(CCCurDataSSplitO,2),EvenOddMap(1,:,:,:,:,:,:,:,:,:,:));
        Cat=CombineDataDim(CCCurDataSSplitE,ifft1cg(CCCurDataSCOC,2));
        nDims=length(size(squeeze(Cat)));
        CmplxDataEO=CombineDims(squeeze(Cat),[1 nDims]);
        
        DataFinal=fft1cg(CmplxDataEO,2);
        
        SensFull=MultMatTensor(FinalAFull',DataFinal);
        FinalB0Recon(:,:,:,SliI)=SensFull;
    end
end