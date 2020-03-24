function ResultRate=Sim_Ephotosynthesis
InputTable=importdata('InputEvn.txt');
EvnInput=InputTable.data;
Cai=EvnInput(1);
Lii=EvnInput(2);
SucPath=EvnInput(3);%if SucPath=1 with sucrose synthesis(original)///if SucPath=0 No Sucrose synthesis,only T3P output
ATPCostTable=importdata('InputATPCost.txt');
ATPCost=ATPCostTable.data;
ResultRate=trDynaPS_Drive(Cai,Lii,ATPCost,SucPath,1, 1);
dlmwrite('OutputRate.txt',ResultRate)
end

