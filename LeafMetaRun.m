

global pcfactor;  
ProteinTotalRatio=0.973;
pcfactor=1/ProteinTotalRatio;
inE=importdata('MeM_input5_0.txt');
Einput=inE.data;
Edata=importdata('Einput7.txt');
Eio=Edata.data(:,1);
MetaOnly=0;% if MetaOnly=1 run Metabolic model

WeatherTemp=25;
Air_CO2=400;
WeatherRH=0.6;
WeatherWind=5;
Radiation_PAR=470*0.85*0.85;%%%%
Convert=1E6/(2.35E5); %Convert W m^{-2} to u moles m^{-2} s^{-1}
Radiation_NIR=0;
Radiation_LW=0;
PhotosynthesisType=1.1;
Vcmax25=100;
Jmax25=200;
GRNC=1;

if MetaOnly==1
CO2i=Air_CO2*0.7; % intercellular CO2 
PPFDi=Radiation_PAR*Convert;
NetAssimilation=EPS_Drive_GRNs(Einput,CO2i,PPFDi,WeatherTemp,GRNC,0,Eio);
else
LeafResult=Leaf(WeatherRH,WeatherTemp,Air_CO2,WeatherWind,Radiation_PAR,Radiation_NIR,Radiation_LW,PhotosynthesisType,Vcmax25,Jmax25,GRNC,Einput,Eio);
Ci=LeafResult(1);
NetAssimilation=LeafResult(2);
Gs=LeafResult(3);
LeafTemperature=LeafResult(4);
Transpirationi=LeafResult(5);
end
fileID = fopen('LeafmetaOut.txt','w');
fprintf(fileID,'%6s\n','A');
fprintf(fileID,'%6.2f\n',NetAssimilation);
fclose(fileID);