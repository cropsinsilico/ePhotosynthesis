function CM_Drive_GRN(inputfile1,inputfile3)
%CM_Drive_GRN('MeM_input1_2.txt','MeM_input3.txt' )
global GRNC;
global GRNT;
indata1= importdata(inputfile1);
indata3= importdata(inputfile3);
GRN_data=indata1.data; 
Env_data=indata3.data; 

CO2in = Env_data(1);
Liin = Env_data(2);
Tpin = Env_data(3);
GRNC = Env_data(4);
GRNT = Env_data(5);

global VfactorC;
global VfactorT;

VfactorC=GRN_data;

Arate=MeM_run(Liin,CO2in,Tpin)
end