InputTable=importdata('InputEvn.txt');
EvnInput=InputTable.data;
Cai=EvnInput(1);
SucPath=EvnInput(3);%if SucPath=1 with sucrose synthesis(original)///if SucPath=0 No Sucrose synthesis,only T3P output
ATPCostTable=importdata('InputATPCost.txt');
ATPCost=ATPCostTable.data;
Begin = 1;
fin = SYSInitial(Begin);
global tglobal;
time=tglobal;%simulation time
CA=[100,150,200,250,300,400,500,600,800,1200,1500,1800];
global VmaxAdj;
VmaxAdj=1.36;%adjust enzyme activity

for i=1:12  
    i
Lii=1500;
Cai=CA(i);
ResultRate(:,i)=trDynaPS_Drive(Cai,Lii,ATPCost,SucPath,1, 1,time);

end
ACI(:,2)=0.7*CA';
ACI(:,1)=(ResultRate(1,:)-0.5*ResultRate(2,:))';