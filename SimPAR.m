% PCfactorMri=[0.964,1.002,0.960,0.984,0.977,1.003,0.981,1.007,0.986,1.011,0.998,1.021,1.012,1.033,0.972,0.999];
% inE=importdata('MeM_input8.txt');
% PCfactorMri=[0.700,0.687,1.01,1.037,0.961,0.980];
% inE=importdata('MeM_input9.txt');
% PCfactorMri=[0.897,0.9225,0.821,0.826,0.704,0.695];
% inE=importdata('MeM_input10.txt');
% 
% PCfactorMri=[0.973,1];

global pcfactor;  
ProteinTotalRatio=0.973;
pcfactor=1/ProteinTotalRatio;
inE=importdata('MeM_input5_0.txt');
Einput=inE.data;
Edata=importdata('Einput7.txt');
Eio=Edata.data(:,1);


for i=1:25
WeatherTemp=25;
Air_CO2=400;
WeatherRH=0.6;
WeatherWind=5;
Radiation_PAR=20*i-10;
Radiation_PARi(i)=Radiation_PAR;
Radiation_PARa=Radiation_PAR*0.85*0.85;
Radiation_NIR=0;
Radiation_LW=0;
PhotosynthesisType=1.1;
Vcmax25=100;
Jmax25=200;
GRNC=0;
LeafResult=Leaf(WeatherRH,WeatherTemp,Air_CO2,WeatherWind,Radiation_PARa,Radiation_NIR,Radiation_LW,PhotosynthesisType,Vcmax25,Jmax25,GRNC,Einput,Eio);
%Result=Leaf(0.6,28,5,400,0,0,1)

Cii(i)=LeafResult(1);
NetAssimilationi(i)=LeafResult(2);
Gsi(i)=LeafResult(3);
LeafTemperaturei(i)=LeafResult(4);
Transpirationi(i)=LeafResult(5);
end
Radiation_PAR=0;

for j=1:25
WeatherTemp=25;
Air_CO2=400;
WeatherRH=0.6;
WeatherWind=5;
Radiation_PAR=20*j-10;
Radiation_PARa=Radiation_PAR*0.85*0.85;
%Radiation_PARi(i)=Radiation_PAR;
Radiation_NIR=0;
Radiation_LW=0;
PhotosynthesisType=1.1;
Vcmax25=100;
Jmax25=200;
GRNC=1;
LeafResult2=Leaf(WeatherRH,WeatherTemp,Air_CO2,WeatherWind,Radiation_PARa,Radiation_NIR,Radiation_LW,PhotosynthesisType,Vcmax25,Jmax25,GRNC,Einput,Eio);
%Result=Leaf(0.6,28,5,400,0,0,1)

HCii(j)=LeafResult2(1);
HNetAssimilationi(j)=LeafResult2(2);
HGsi(j)=LeafResult2(3);
HLeafTemperaturei(j)=LeafResult2(4);
HTranspirationi(j)=LeafResult2(5);
end

Radiation_PARi=Radiation_PARi*1E6/(2.35E5); %Convert W m^{-2} to u moles m^{-2} s^{-1}
figure;
subplot(2,2,1); plot(Radiation_PARi,NetAssimilationi,'r');hold on; plot(Radiation_PARi,HNetAssimilationi,'k');xlabel('PAR');ylabel('Net CO2 assimilation');
subplot(2,2,2); plot(Radiation_PARi,Transpirationi,'r');hold on; plot(Radiation_PARi,HTranspirationi,'k');xlabel('PAR');ylabel('Transpiration');
subplot(2,2,3); plot(Radiation_PARi,Gsi,'r');hold on; plot(Radiation_PARi,HGsi,'k'); xlabel('PAR');ylabel('Gs');
subplot(2,2,4); plot(Radiation_PARi,LeafTemperaturei,'r');hold on; plot(Radiation_PARi,HLeafTemperaturei,'k');xlabel('PAR');ylabel('Leaf Temperature');