function Out=grms(Data,Dims)
Data=abs(Data).^2;
S=sort(Dims,'descend');
Out=Data;
for i=1:numel(Dims)
    Out=mean(Out,S(i));
end
Out=sqrt(Out);