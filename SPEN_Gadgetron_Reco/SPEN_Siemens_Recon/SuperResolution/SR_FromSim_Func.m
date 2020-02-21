function [FinalA_Sim,B0Extend]=SR_FromSim_Func(spinsA, FOV, Rval, Nacq, EchoSpacing,NpeTrg,BWEnc, Te, PTa2Delay,XX,PChirp,PPERephaseBefore,PPERephaseAfter,PHalfBlip,PBlip,LPE)
PEShiftx=0;
R=[spinsA.r];
B0Norm=zeros(NpeTrg,1); % SpenSimulation_InPlaneClean; %SpenSim_SR_function
B0Extend=interp1(linspace(-FOV(2)/2,FOV(2)/2,NpeTrg),B0Norm',R(2,:)+PEShiftx,'linear','extrap');

% Using Data
N=size(B0Extend);
SpinsPerVox=N(2)/NpeTrg;
R=[spinsA.r];

t0_fraction = (R(2,:))/FOV(2) + B0Extend/BWEnc; % [-0.5 0.5]
t0_01=t0_fraction+0.5; % [0 1]
t0_us=t0_01*Te;
t0_us_centered=t0_fraction*Te;
% B0Points=B0Extend;%0.M0X1*(B0Extend-1)/2;
for i=1:N(2)
    spinsA(i).B0=B0Extend(i);
end

MStart=[spinsA.M];

PhiStart=atan2(MStart(2,:),MStart(1,:));

M0X=B0Extend*0+1;

PhiAfterDelay=PhiStart-2*pi*PTa2Delay.tp*(B0Extend+R(2,:)*PTa2Delay.Gy(1));

PhiAfterPrephase=PhiAfterDelay-2*pi*PPERephaseBefore.tp*(B0Extend+R(2,:)*PPERephaseBefore.Gy);

ChirpPhiAt_t0= -2*pi*(t0_us_centered.^2) * Rval/(Te.^2);
ChirpPhiAt_t0=ChirpPhiAt_t0*1.5;

PhiAtFlipTime=PhiAfterPrephase-2*pi*(B0Extend+R(2,:)*PChirp.Gy(1)).*t0_us/1000;
PhiJustAfterFlip=2*ChirpPhiAt_t0-PhiAtFlipTime;
PhiAfterChirp=PhiJustAfterFlip-2*pi*(B0Extend+R(2,:)*PChirp.Gy(1)).*(1-t0_01)*Te/1000;

PhiX=PhiAfterChirp;%+ChirpPhiAt_t0;

PhiX=PhiX-2*pi*PPERephaseAfter.tp*(B0Extend+R(2,:)*PPERephaseAfter.Gy);
PhiX=PhiX-XX*2*pi*PHalfBlip.tp*(B0Extend+R(2,:)*PHalfBlip.Gy);

BlipMoment=PBlip.tp*PBlip.Gy;
SigY=NaN(NpeTrg,Nacq);
%
M0XReshap=reshape(M0X,[SpinsPerVox,NpeTrg]);
for i=1:Nacq
    PhiXReshap=reshape(PhiX,[SpinsPerVox,NpeTrg]);

    SigY(:,i)=mean(M0XReshap.*exp(1i*PhiXReshap));
    PhiX=PhiX-2*pi*EchoSpacing*B0Extend-2*pi*R(2,:)*BlipMoment;
end
% SR=mean(SigY,1);

FinalA_Sim=flip(SigY')/Nacq *LPE;


