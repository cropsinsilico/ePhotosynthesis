
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
global VmaxAdj;
VmaxAdj=1.36;%adjust enzyme activity

for i=1:30
i
Lii=50*i;

ResultRate=trDynaPS_Drive(Cai,Lii,ATPCost,SucPath,1, 1,time);

if ResultRate(9)>5e-05
     warning('ODE model does not reach a steady state: The 1st time ');
     ResultRate=trDynaPS_Drive(Cai,Lii,ATPCost,SucPath,1, 1,time*10);% if the ODE model can't reach steady state, increase the simulation time.
     if ResultRate(9)>5e-05
         warning('ODE model does not reach a steady state: The 2nd time (the last time) ');
         ResultRate=0;
     end

end
ResultRatei(:,i)=ResultRate;
end