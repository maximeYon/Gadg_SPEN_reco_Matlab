function FinalA_R4=SuperResolution_VariableSize_Sim_V3(Nacq,ParamSR,R,WeightingFactor,nShots,ChirpAmpExtraFactor)
% Work for different R value
% T2S_Estimate_ms=50; in brain
% Chirp use to not have a big influence on the quality of the
% reconstruction
% FOV=[0 mrprot.sSliceArray.asSlice{1}.dPhaseFOV]; % mm

LPE=ParamSR.LPE;
PEShift=ParamSR.PEShift;
rvalue=ParamSR.rvalue;
EchoSpacing=ParamSR.EchoSpacing;
if exist('ParamSR.Shift_vx')
    Shift_vx=ParamSR.Shift_vx;
else
    Shift_vx=0;
end

Scl=1000;
if(nargout>1)
    Scl=100;
end

% Shift_vx=0; 

FOV=[0 LPE*10]; % mm
EchoSpacing_ms=EchoSpacing*0.001;
DwellTime=10; % us
if(R>9)
    Nvox=R;
else
    Nvox= Nacq*R;
end
M0=ones(1,Nvox*Scl);
AcqBWToEncBWRatio=1;
cs_KHz=0;

Ta=EchoSpacing*Nacq/nShots; % ms
Te=Ta/2; % that'll change for non-fully-refocused

PEShiftRelativeToFOV=-PEShift/FOV(2);
PEShiftx=FOV(2)*PEShiftRelativeToFOV; % mm


[spinsA, PTa2Delay,PChirp,PPERephaseBefore,PHalfBlip,PPERephaseAfter,PBlip]=... % pulse_preparation_sim %%Modify not tested. SuperResolution_VariableSize_Sim_V3
    Pulse_Prep_Function(M0,Nacq,Nvox,EchoSpacing_ms,Ta, Te, rvalue, FOV, AcqBWToEncBWRatio,-PEShiftx,0); 

% Add phase to chirp according to PEShift
nSteps=Te/DwellTime;
BW_kHz=1000*rvalue/Te;
FreqShift_kHz=BW_kHz*PEShift/FOV(2);
PointTime_us=((1:nSteps)-nSteps/2)*DwellTime;
AddedShiftPhase=(FreqShift_kHz/1000)*PointTime_us;
PChirp.RFphase=PChirp.RFphase+AddedShiftPhase*2*pi;

% Correction factor for the Chirp
PChirp.RFamp=PChirp.RFamp*ChirpAmpExtraFactor;

% Another correction
BWEnc=rvalue*1000/Te; % kHzSuperResolution_VariableSize_Sim_V3
PEShiftRelativeToFOV=-PEShift/FOV(2);
BlipCorrection=2-2*Nacq*PEShiftRelativeToFOV;
% Final A
FinalA_R4=SR_FromSim_Func(spinsA, FOV, rvalue, Nacq, EchoSpacing_ms,Nvox,BWEnc, Te, PTa2Delay,BlipCorrection,PChirp,PPERephaseBefore,PPERephaseAfter,PHalfBlip,PBlip,LPE);

% Limitation for high frequency
[X, Y]=meshgrid(linspace(1,Nacq,size(FinalA_R4,2)),linspace(1,Nacq,size(FinalA_R4,1)));
XY=abs(X-Y);

EchoSpacing_us=EchoSpacing_ms*1000;
FinalA_R4=FinalA_R4.*exp(-XY*((EchoSpacing_us/1e3)/nShots)/WeightingFactor);

