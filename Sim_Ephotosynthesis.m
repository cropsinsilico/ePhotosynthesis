function ResultRate=Sim_Ephotosynthesis
InputTable=importdata('InputEvn.txt');
EvnInput=InputTable.data;
Cai=EvnInput(1);
Lii=EvnInput(2);
SucPath=EvnInput(3);%if SucPath=1 with sucrose synthesis(original)///if SucPath=0 No Sucrose synthesis,only T3P output
ATPCostTable=importdata('InputATPCost.txt');
ATPCost=ATPCostTable.data;
NADPHCostTable=importdata('InputNADPHCost.txt');
NADPHCost=NADPHCostTable.data;
Begin = 1;
fin = SYSInitial(Begin);
global tglobal;
time=tglobal;%simulation time
global VmaxAdj;
VmaxAdj=1.36;%adjust enzyme activity


ResultRate=trDynaPS_Drive(Cai,Lii,ATPCost,NADPHCost,SucPath,1, 1,time);
ResultRate(9)
if ResultRate(9)>5e-05
     warning('ODE model does not reach a steady state: The 1st time ');
     ResultRate=trDynaPS_Drive(Cai,Lii,ATPCost,NADPHCost,SucPath,1, 1,time*10);% if the ODE model can't reach steady state, increase the simulation time.
     if ResultRate(9)>5e-05
         warning('ODE model does not reach a steady state: The 2nd time (the last time) ');
         ResultRate(:,1)=0;
     end

end

Output(1)="Light intensity,Vc,Vo,VPGA,VT3P,Vstarch,Vt_glycerate,Vt_glycolate,A_slope,Max_Meta_slope";
Output(2)=Lii+","+ResultRate(1)+","+ResultRate(2)+","+ResultRate(3)+","+ResultRate(4)+","+ResultRate(5)+","+ResultRate(6)+","+ResultRate(7)+","+ResultRate(8)+","+ResultRate(9);
outfile = fopen("OutputRate.txt","w");
fprintf(outfile,Output(1));
fprintf(outfile,"\n");
fprintf(outfile,Output(2));
fclose(outfile);
%%%output flux file
global FluxTR;
FluxName=importdata('FluxName.txt');
vect_out(2,:) = num2cell(FluxTR);
vect_out(1,:) = FluxName;
fid = fopen('OutputFluxT.txt','wt');
fprintf(fid,'%-14s  %d\n',vect_out{:});
fclose(fid);
end

