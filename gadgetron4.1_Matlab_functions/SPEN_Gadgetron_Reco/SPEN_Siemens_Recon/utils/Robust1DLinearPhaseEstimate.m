function [Fld, BestX, BestCost]=Robust1DLinearPhaseEstimate(In,dim,DoOptimization)
if(nargin==2)
    DoOptimization=true;
end
OtherDims=setdiff(1:ndims(In),dim);
In=permute(In,[dim OtherDims]);
In=double(In);
W=abs(In);
WForD=(abs(In(1:end-1,:))+abs(In(2:end,:)))/2;
OnlyP=exp(1i*angle(In));
WD=WForD.*OnlyP(2:end,:)./OnlyP(1:end-1,:);
DPEst=angle(sum(WD(:)));
Sz=size(In);
OnlyLinearFldCol=((1:Sz(1))*DPEst).';
OnlyLinearFld=W.*exp(1i*repmat(OnlyLinearFldCol,[1 Sz(2:end)]));
% est bias
Bias=OnlyLinearFld(:).'/max(In(:).',eps);
LinearFldWBiasCol=OnlyLinearFldCol-angle(Bias);
LinearFldWBias=exp(1i*repmat(LinearFldWBiasCol,[1 Sz(2:end)]));

% optimize
x0=[angle(Bias) DPEst];
ColFunc=@(x) ((1:Sz(1))*x(2)).'-x(1);
FldFunc=@(x) W.*exp(1i*repmat(ColFunc(x),[1 Sz(2:end)]));
CostFunc=@(x) gsum(abs(FldFunc(x)-In));
% [BestX BestCost]=fminsearch(CostFunc,[0 0])
BestX=x0;
BestCost=CostFunc(x0);
Fld=FldFunc(x0);
if(~DoOptimization)
    return;
end
[BestX, BestCost]=fminsearch(CostFunc,x0);
Fld=FldFunc(BestX);
Fld=ipermute(Fld,[dim OtherDims]);