function out = PEShiftingDataset(Data,TotalPEShift,LPE,nShots) 
%input data in fourier domaine

CurData=fft1cg(Data,2);
tmp=PartitionDimInterleaved(CurData,1,nShots);

NpeS=size(tmp,1);
%         tmp=RepDotMult(tmp,         exp(1i*(2*pi/LPE * (TotalPEShift/10)/nShots) * (1:NpeS) ).');
tmp=RepDotMult(tmp,         exp(1i*(2*pi/LPE * (TotalPEShift/10)) * (1:NpeS) ).');

if nShots==1
    CurDataShift=tmp;
else
    CurDataShift=CombineDims(tmp,[1 length(size(tmp))]);
end
out=ifft1cg(CurDataShift,2);