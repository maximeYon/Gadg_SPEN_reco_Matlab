function Out=gsum(in,dims)
if(~exist('dims','var'))
    Out=sum(in(:));
    return;
end
Out=sum(in,dims(1));
for i=2:length(dims)
    Out=sum(Out,dims(i));
end