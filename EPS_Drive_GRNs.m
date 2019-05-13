function Arate=EPS_Drive_GRNs(input1,CO2i,PPFDi,Tempi,Gc,GT,Einput)
global GRNC;
global GRNT;
global pcfactor;
GRN_data=input1*pcfactor;

CO2in = CO2i;
Liin = PPFDi;
Tpin = Tempi;
GRNC = Gc;
GRNT = GT;

global EnzymeAct;
EnzymeAct=Einput(1:27)/30;%unit change
global Jmax;
Jmax=EnzymeAct(27);
global BFVmax;
BFVmax=Einput(28:45);
global FIVmax;
FIVmax=Einput(46:66);

global cATPsyn;
global CPSi;
global cNADPHsyn;
global cpsii;
if GRNC==1
cATPsyn=GRN_data(34);%1.0447;%1.01866 WY201803
CPSi=GRN_data(35);%1.0131;% 1.0237 WY201803
cNADPHsyn=GRN_data(37);%1.094468408;%1.0388 WY201803
cpsii=GRN_data(36);%1.0169;% 1.0129;%WY201803
end
if GRNC==0
cATPsyn=1;%1.01866 WY201803
CPSi=1;% 1.0237 WY201803
cNADPHsyn=1;%1.0388 WY201803
cpsii=1;% 1.0129;%WY201803   
end
global VfactorC;
global VfactorT;

VfactorC=GRN_data(1:33);

Arate=EPS_Drive(Liin,CO2in,Tpin);
end