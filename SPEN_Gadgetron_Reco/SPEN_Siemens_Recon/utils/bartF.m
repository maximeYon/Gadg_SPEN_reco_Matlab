function FinalRecon=bartF(DataC,FinalAC,CombinedSensPhaseDiff,ReconIndx,TVW)

DataP=permute(DataC,[4 5 2 3 6 1]);
%             DataP=DataP+randn(size(DataP))*grmss(DataP)*0.1;
            
SRP=repmat(permute(FinalAC,[5 2 4 3 6 1]),[1 1 size(DataC,2) 1 1 1]);
SensP=permute(CombinedSensPhaseDiff,[4 1 2 3]);
FinalRecon=squeeze(bart(ReconIndx,DataP,SensP,SRP*TVW));