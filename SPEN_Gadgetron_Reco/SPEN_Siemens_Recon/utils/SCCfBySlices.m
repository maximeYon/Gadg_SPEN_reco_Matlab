function [sccmtx, Energy]=SCCfBySlices(Im,ncalib)
if(nargin<2)
    ncalib = 20; % use 24 calibration lines to compute compression
end
nSlices=size(Im,4);
for s=1:nSlices
    DATA=fft2cg(Im(:,:,:,s));
    [sx,sy,Nc] = size(DATA);
    calib = crop(DATA,[ncalib,sy,Nc]);
    sccmtx(:,:,s) = calcSCCMtx(calib);
    if(nargout>1)
        SCCDATA = CC(DATA,sccmtx(:,:,s));
        tmp=grmss(SCCDATA,1:2);
        Energy(:,s)=tmp./sum(tmp);
    end
end