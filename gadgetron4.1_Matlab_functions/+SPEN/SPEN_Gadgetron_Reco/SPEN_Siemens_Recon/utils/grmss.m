function Out=grmss(Data,Dims)
if(nargin<2)
    Dims=1:ndims(Data);
end
Out=squeeze(grms(Data,Dims));