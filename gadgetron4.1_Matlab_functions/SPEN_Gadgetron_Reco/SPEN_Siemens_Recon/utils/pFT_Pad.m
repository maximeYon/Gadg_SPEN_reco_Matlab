function [padded_CmplxDataSh]=pFT_Pad(CmplxDataSh,pFT)
nRO=size(CmplxDataSh,2);
NeededToPad=ceil(nRO/pFT-nRO);
%padded_CmplxDataSh=[CmplxDataSh zeros([size(CmplxDataSh,1) NeededToPad gsize(CmplxDataSh,3:8)])];
padded_CmplxDataSh=[zeros([size(CmplxDataSh,1) NeededToPad  gsize(CmplxDataSh,3:7)]) CmplxDataSh]; % for 3D SPEN pFT is in other direction
end


% test=grms(CmplxData,[1 3 4 5 6]);
% figure,plot(test)
