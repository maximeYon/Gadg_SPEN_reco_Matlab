function [FdataCompressed, gccmtx_aligned_final, Energy]=CoilCompression(FDataToCompress,ncc)
% Details :
% Perform coil compression based on work of Zhang et. al MRM 2013;69(2):571-82.
% Input (FDataToCompress) should be raw data with shots already combined.
% ncc: target number of channels
nChannels=size(FDataToCompress,3);

if nChannels>ncc
    disp('SCCing');
    [gccmtx_aligned_final, Energy]=SCCfBySlices(FDataToCompress(:,:,:,:,1,1,1,1));
    FdataCompressed=ApplySCCBySlices(FDataToCompress,gccmtx_aligned_final,ncc);
else
    disp('Not appliting CC, same ncc');
    FdataCompressed=FDataToCompress;
    gccmtx_aligned_final=eye(ncc);
    Energy=zeros(1,ncc);
end