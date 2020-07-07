function  res = SPEN_Sens_SR4_MB(PEFactor,Sens,SRMat,mode)

%res = p2DFT(mask,imSize [ ,phase,mode])
%
%
%	Implementation of partial Fourier operator.
%	
%	input:
%			mask - 2D matrix with 1 in entries to compute the FT and 0 in ones tha
%				are not computed.
%			imSize - the image size (1x2)
%			Sens - Multi slice sensitivity maps [PE, RO, Slice]
%			mode - 1- real, 2-cmplx
%
%	Output:
%			The operator
%
%	(c) Michael Lustig 2007

if nargin <4
	mode = 2; % 0 - positive, 1- real, 3-cmplx
end

% res.PEFactor=PEFactor;
imSize=gsize(Sens,1:3);
res.adjoint = 0;
res.imSize = imSize;
res.dataSize = imSize;
res.dataSize(1)=round(res.dataSize(1)/PEFactor);
if(numel(mode)==1)
    mode(2)=0;
end
res.DoROFT(1)=mode(2);
if(res.DoROFT(1))
    res.DoROFT(2)=mode(3);
end
if(numel(res.dataSize)<3)
    res.dataSize(3)=1;
end
res.Sens = Sens;
res.mode = mode(1);
res.nChannels=size(Sens,3);
res.nBands=size(Sens,4);
res.invSens=conj(permute(Sens,[1 2 4 3]));
res.SRMat=SRMat;
res = class(res,'SPEN_Sens_SR4_MB');