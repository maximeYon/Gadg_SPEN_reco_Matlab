function WURSTEnvelope=CreateWURSTEnvelope(n,WURSTn)

WURSTEnvelope=1-abs(sin(linspace(-pi/2,pi/2,n)).^WURSTn);