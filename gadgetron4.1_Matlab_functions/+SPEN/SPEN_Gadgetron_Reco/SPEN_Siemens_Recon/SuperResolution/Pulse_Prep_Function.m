function [spinsA, PTa2Delay,PChirp,PPERephaseBefore,PHalfBlip,PPERephaseAfter,PBlip]=Pulse_Prep_Function(M0,Nacq,Nvox,EchoSpacing,Ta_ms, Te, Rval, FOV,AcqBWToEncBWRatio,PEShiftSim,cs_KHz)
% Inputs

HalfTa=Ta_ms/2;
DwellTime=10;
N=size(M0);
SpinsPerVox=N(2)/Nvox;
%% OK, preparing the "Pulses"

%Delay
PTa2Delay=PulseCreateHard(HalfTa, 0, 0, 100) ;

%Chirp
PChirp=CreateChirpPulseG(Te,Rval,FOV,DwellTime); % 90 chirp
PChirp.RFamp=PChirp.RFamp*2.25; % Becomes 180 chirp
PChirp.RFamp=PChirp.RFamp.*CreateWURSTEnvelope(numel(PChirp.RFamp),40);

PrePhaserMoment=AcqBWToEncBWRatio*PChirp.tp*PChirp.Gy(1);
HalfPrePhaserMoment=PrePhaserMoment/2;

PPERephaseBefore.tp=0.01; 
PPERephaseBefore.RFamp=0; PPERephaseBefore.RFphase=0; PPERephaseBefore.Gx=0; PPERephaseBefore.Gz=0;
PPERephaseBefore.Gy=HalfPrePhaserMoment/PPERephaseBefore.tp; % PChirp.Gy(1)*100*2/(Nacq);

PPERephaseAfter.tp=0.01; 
PPERephaseAfter.RFamp=0; PPERephaseAfter.RFphase=0; PPERephaseAfter.Gx=0; PPERephaseAfter.Gz=0;
PPERephaseAfter.Gy=-HalfPrePhaserMoment/PPERephaseAfter.tp; % PChirp.Gy(1)*100*2/(Nacq);

TotalBlipMoment=AcqBWToEncBWRatio*PChirp.tp*PChirp.Gy(1)*2/(Nacq);
PBlip.tp=0.01; % 10 us; % Should be replaced with some minimal value, e.g. 10us
PBlip.RFamp=0; PBlip.RFphase=0; PBlip.Gx=0; PBlip.Gz=0;



PBlip.Gy=TotalBlipMoment/PBlip.tp; % PChirp.Gy(1)*100*2/(Nacq);
PHalfBlip=PBlip;
PHalfBlip.Gy=PHalfBlip.Gy/2;

PRONoGradForEchoSpacing=PBlip;
PRONoGradForEchoSpacing.Gy=0;
PRONoGradForEchoSpacing.tp=EchoSpacing-PBlip.tp; % e.g. for echospacing of 700us 

%Spin preparation

spinsA=InitSpinsRelax(0, N(2), 1, [0; 0; 1], Inf, Inf, 1);

BinEdges=linspace(-FOV(2)/2,FOV(2)/2,Nvox+1);
BinCenters=(BinEdges(1:end-1)+BinEdges(2:end) )/2;
BinSize=BinEdges(2)-BinEdges(1);

AllLocs=repmat(BinCenters,[SpinsPerVox 1]);
AllLocs=AllLocs+(rand(SpinsPerVox,Nvox)-0.5)*BinSize;

tmp=linspace(-FOV(2)/2,FOV(2)/2,N(2)+1);
AllLocs(:)=tmp(1:end-1);
AllLocs=AllLocs+(rand(SpinsPerVox,Nvox)-0.5)*BinSize/SpinsPerVox;


for i=1:N(2)
    spinsA(i).r(1)=0; % mm
%     yN=   (       rand-0.5     +    i-(   N(2)/2     )   )    /   N(2)              ; % a.u. - between -0.5 to 0.5 along the FOV
%     spinsA(i).r(2)=yN*FOV(2)-PEShiftSim; % mm
    spinsA(i).r(2)=AllLocs(i)-PEShiftSim; % mm
    spinsA(i).r(3)=0; % mm
    spinsA(i).M0=M0(i); % 1; % a.u.
    spinsA(i).M=[spinsA(i).M0; 0; 0]; % a.u.  % we start with magnetization on x (on the plane)
    spinsA(i).cs=cs_KHz; % 1; % kHz
end






%% TRASH 
% 1.
% 
% t=linspace(-Te/2,Te/2,nSteps);
% phi= 2*pi*(t.^2) * Rval/(Te.^2);

% nSteps=Te/DwellTime;

% MaxPha = pi*Rval/2;
% a_rad2mm2=  MaxPha/(FOV(2)/2);%2*pi*Rval/FOV(2); %
% BW_sp = 2*a_rad2mm2;
% PredictedRes_mm = 4/BW_sp; % = 2*FOV/pi*Rval
% MaxNpe=FOV(2)/PredictedRes_mm;
% 

% MaxNpe=FOV(2)*BW_sp/4 = FOV(2)*a_rad2mm2/2 = FOV(2)*(pi*Rval/FOV(2))/2 = pi*Rval/2;