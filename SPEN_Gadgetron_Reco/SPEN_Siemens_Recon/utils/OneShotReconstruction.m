% This reconstruction is for one shot data
%
% The workflow is:
%
% 1) Compress data to one and nccx channes;
% 2) Select one slice, one bval;
% 3) Split data and SR matrix into even and odd data subsets;
% 4) Reconstruct even and odd images from data compressed to one channel
%      with the resolution of final image, determine even-odd correction from
%      phase difference between even and odd data using 1D phase fitting;
% 5) Split data compressed to nccx channels into even and odd data sets,
% apply even-odd correction determined for data compressed to one channel
% 6) Determine sensitivity maps from even-odd corrected data and
% reconstructed by simple multiplication with Superresolution matrix.
% 7) Reconstruct image with some L2 regularization.
%% Coil Compressing
CCCurCmplxDataSx=ApplySCCBySlices(CmplxDataSh(:,:,:,:,:,:,:),gccmtx_aligned_finalx,nccx); % Compressing data to nccx channels
CCCurCmplxDataS1=ApplySCCBySlices(CmplxDataSh(:,:,:,:,:,:,:),gccmtx_aligned_finalx,1); % Compressing data to one channel

clear SensFull
clear FinalB0Recon
clear FinalRep


for bval=1:nWeight
    for SliI=1:nSlices
        %% Even Odd
        CurDataEO=CCCurCmplxDataS1(:,:,:,SliI,:,:,bval);
        
        CCCurDataSSplitE=CurDataEO(1:2:end,:,:,:,:,:,:,:); % Select even data
        CCCurDataSSplitO=CurDataEO(2:2:end,:,:,:,:,:,:,:); % Select odd data
        
        FinalATrySplitE=FinalAFull(1:2:end,:,:,:,:); % Select even lines from SR matrix
        FinalATrySplitO=FinalAFull(2:2:end,:,:,:,:); % Select even lines from SR matrix
        SensShotCorrection=ones([size(FinalATrySplitE,2) gsize(CCCurDataSSplitE,[2 3 4 5 6 7 8 9 10])]);
        
        FinalE=RunfnlCgIterationsf(1000+1,FinalATrySplitE,11,squeeze(SensShotCorrection(:,:,:,:,1,1,:)),TVW_L2.EO,squeeze(fft1cg(CCCurDataSSplitE(:,:,:,:,:,1,:),2)),2,struct('ShowFig',false)); % Image Even
        FinalO=RunfnlCgIterationsf(1000+1,FinalATrySplitO,11,squeeze(SensShotCorrection(:,:,:,:,1,1,:)),TVW_L2.EO,squeeze(fft1cg(CCCurDataSSplitO(:,:,:,:,:,1,:),2)),2,struct('ShowFig',false)); % Image Odd
        %%
        CurW=abs(FinalE);
        SRPhaseDiff=angle(FinalE)-angle(FinalO); % phase difference between odd and even
        EOWithW=CurW.*exp(1i*SRPhaseDiff); % some weight on the phase
        [ROEvenOddCorrection_bvalu, BestX, BestCost] =Robust1DLinearPhaseEstimate(EOWithW,2); % do a 1D linear phase fitting on the phase that will later be applied to all readout lines
        EvenOddMap=exp(+1i*angle(repmat(ROEvenOddCorrection_bvalu(1,:,:,:,:,:,:,:),[Npe 1 1 1 1 ])));
        
         %% Apply EO to data compressed to 1 channel 
        CCCurDataSCOC_1ch=RepDotMult( fft1cg(CCCurDataSSplitO,2),EvenOddMap(1,:,:,:,:,:,:,:,:,:,:)); % Apply even-odd correction
        CmplxDataEO_1ch=CombineDataDim(CCCurDataSSplitE,ifft1cg(CCCurDataSCOC_1ch,2)); % Recombine the data
        DataFinal_1ch=fft1cg(CmplxDataEO_1ch,2); % Do FT along RO dimension

        
        %% Apply EO to multicoil data
        CCCurDataSSplitE=CCCurCmplxDataSx(1:2:end,:,:,SliI,:,:,:,:); % Select even data from data split to nccx channels
        CCCurDataSSplitO=CCCurCmplxDataSx(2:2:end,:,:,SliI,:,:,:,:); % Select even data from data split to nccx channels
        
        CCCurDataSCOC=RepDotMult( fft1cg(CCCurDataSSplitO,2),EvenOddMap(1,:,:,:,:,:,:,:,:,:,:)); % Apply even-odd correction
        CmplxDataEO=CombineDataDim(CCCurDataSSplitE,ifft1cg(CCCurDataSCOC,2)); % Recombine the data
        
        DataFinal=fft1cg(CmplxDataEO,2); % Do FT along RO dimension
        
        for NCha=1:nccx
            SensFull(:,:,NCha)=MultMatTensor(FinalAFull',DataFinal(:,:,NCha,:,1,1,1)); % Get image by multiplication with SR matrix
        end
        
        % Determine Sensitivity maps
        if nccx==1
            Sens=ones(size(SensFull));
        else
            Sens=RunESPIRiTForSensMaps( SensFull,24,[6 6]);
        end
        
        % Reconstruct image
        for RIndx=1:nReps
            L2image(:,:,SliI,RIndx,bval,:)=RunfnlCgIterationsf(1000+1,FinalAFull,2,Sens,TVW_L2.Final,DataFinal(:,:,:,:,:,RIndx,bval),2,struct('ShowFig',false));
        end
                
    end
end