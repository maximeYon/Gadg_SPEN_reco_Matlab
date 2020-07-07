function Out=ExtendStruct(Src, Full)
% function Out=ExtendStruct(Src, Full)
F=fieldnames(Src);
ToAdd=union(setdiff(fieldnames(Full)',fieldnames(Src)'),F(gIsEmpty(struct2cell(Src))));
ToAdd=ToAdd(:)';
Out=Src;
for i=ToAdd
    Out.(i{1})=Full.(i{1});
end