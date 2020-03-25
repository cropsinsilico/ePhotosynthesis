function ResultRate=Sim_Ephotosynthesis
InputTable=importdata('InputEvn.txt');
EvnInput=InputTable.data;
Cai=EvnInput(1);
Lii=EvnInput(2);
SucPath=EvnInput(3);%if SucPath=1 with sucrose synthesis(original)///if SucPath=0 No Sucrose synthesis,only T3P output
ATPCostTable=importdata('InputATPCost.txt');
ATPCost=ATPCostTable.data;
ResultRate=trDynaPS_Drive(Cai,Lii,ATPCost,SucPath,1, 1);
Output(1)="Light intensity,Vc,Vo,VPGA,VT3P,Vstarch,Vt_glycerate,Vt_glycolate";
Output(2)=Lii+","+ResultRate(1)+","+ResultRate(2)+","+ResultRate(3)+","+ResultRate(4)+","+ResultRate(5)+","+ResultRate(6)+","+ResultRate(7);
outfile = fopen("OutputRate.txt","w");
fprintf(outfile,Output(1));
fprintf(outfile,"\n");
fprintf(outfile,Output(2));
fclose(outfile);
end

