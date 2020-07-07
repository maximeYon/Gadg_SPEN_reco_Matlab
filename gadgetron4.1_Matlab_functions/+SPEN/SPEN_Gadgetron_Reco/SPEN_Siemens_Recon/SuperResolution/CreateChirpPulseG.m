function Pulse=CreateChirpPulseG(Te,Qval,FOV,DwellTime,Dir)
if(nargin<5)
    Dir=1;
end
nSteps=floor(Te/DwellTime);

BW=1000*Qval/Te; % kHz
R=1000*BW/Te;
Oi=-1*Dir*BW/2; % kHz
Gy=Dir*BW/(FOV(2));
% 2*PI*(Oi * CurTime + (R_chirp/2) * CurTime*CurTime);
TimeN=linspace(0,1,nSteps);
% instFreq=Oi+BW*TimeN
% int(instFreq)=Oi*TimeN+(BW/2)*TimeN^2
% if Freq=Oi, then StartPhi=0, EndPhi=Oi(kHz)*tp(in ms)*2*pi =>
% phi=2*pi*TimeN*Oi*tp=2*pi*tp*int(instFreq)
% Phi=int(instFreq)= Oi(t)+(BW/2)*t^2

Pulse.tp=Te/1000; % ms
Phi=2*pi*Pulse.tp*(Oi*TimeN+(BW/2)*(TimeN.^2));
% Pulse.RFamp=ones(1,nSteps)*(sqrt(Qval)/sqrt(100))*(1/Pulse.tp)*(1/4)*10.5; % kHz  *27.7956
% Pulse.RFamp=ones(1,nSteps)*sqrt(R)*(1/4); % kHz  *27.7956
% sqrt(Qval)/Tp=sqrt(Q)*1000/Te=sqrt(R) !!!!!!
% Pulse.RFamp=ones(1,nSteps)*sqrt(R)*0.27; % kHz  Assaf's paper
Pulse.RFamp=ones(1,nSteps)*sqrt(R)*0.25*1.0598; % kHz  Assaf's paper
Pulse.RFphase=Phi;
Pulse.Gx=zeros(1,nSteps); % kHz/mm
Pulse.Gz=zeros(1,nSteps); % kHz/mm
Pulse.Gy=ones(1,nSteps)*Gy; % kHz/mm