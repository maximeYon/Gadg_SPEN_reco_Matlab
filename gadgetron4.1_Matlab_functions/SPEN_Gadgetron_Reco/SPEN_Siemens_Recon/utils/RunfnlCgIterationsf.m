function im_res=RunfnlCgIterationsf(FigH,SR,nfnlCgIters,SensMFull,TVW,Data,Mode,More)
if(~exist('More','var'))
    More=struct();
end
nChannels=size(SensMFull,3);
nBands=size(SensMFull,4);
ImSize=[gsize(SensMFull,1:2) nBands];
DataSize=size(Data);
PEFactor=ImSize(1)/DataSize(1);
if(size(SR,3)==1)
    SR=repmat(SR,[1 1 nChannels 1]);
end
if(size(SR,4)==1)
    SR=repmat(SR,[1 1 1 nBands]);
end
if(size(SR,1)==ImSize(1))
    SR=SR(1:PEFactor:end,:,:,:);
end
AOdd = SPEN_Sens_SR4_MB(PEFactor,SensMFull,SR,Mode);
OperatorTest(AOdd,ImSize,DataSize,nChannels,1);
param=ExtendStruct(struct('pNorm',2,'TVWeight',TVW,'Itnlim',8,'FT',AOdd,'Verbose',false,'XFM',1,'TV',TVOP_MSlice,'xfmWeight',0),init);
if(exist('More','var'))
    param=ExtendStruct(More,param);
end
param.data =     Data;
res=zeros(ImSize);
if(isfield(param,'WarmStart'))
    res=param.WarmStart;
end
RunFnlViewAmp=1;
if(isfield(More,'RunFnlViewAmp'))
    RunFnlViewAmp=More.RunFnlViewAmp;
end
RunfnlCgIterations;