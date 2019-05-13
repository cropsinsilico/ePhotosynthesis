global pcfactor;  
ProteinTotalRatio=0.973;
pcfactor=1/ProteinTotalRatio;
inE=importdata('MeM_input5_0.txt');
Einput=inE.data;
Edata=importdata('Einput7.txt');
Eio=Edata.data(:,1);
% WeatherTemp=25;
% Air_CO2=400;
% WeatherRH=0.6;
% WeatherWind=5;
% Radiation_PAR=470/2000*2000*0.85;%%%%
% Radiation_NIR=0;
% Radiation_LW=0;
% PhotosynthesisType=1.1;
% Vcmax25=100;
% Jmax25=200;
% GRNC=1;
% LeafResult=Leaf(WeatherRH,WeatherTemp,Air_CO2,WeatherWind,Radiation_PAR,Radiation_NIR,Radiation_LW,PhotosynthesisType,Vcmax25,Jmax25,GRNC,Einput,Eio);
% %Result=Leaf(0.6,28,5,400,0,0,1)
% 
% Cii(i)=LeafResult(1);
% NetAssimilationi(i)=LeafResult(2);
% Gsi(i)=LeafResult(3);
% LeafTemperaturei(i)=LeafResult(4);
% Transpirationi(i)=LeafResult(5);


for i=1:40
WeatherTemp=25;
Air_CO2=10+40*i;
Air_CO2i(i)=Air_CO2;
WeatherRH=0.6;
WeatherWind=5;
Radiation_PAR=470/2000*2000*0.85*0.85;%%%%
Radiation_NIR=0;
Radiation_LW=0;
PhotosynthesisType=1.1;
Vcmax25=100;
Jmax25=200;
GRNC=0;
LeafResult=Leaf(WeatherRH,WeatherTemp,Air_CO2,WeatherWind,Radiation_PAR,Radiation_NIR,Radiation_LW,PhotosynthesisType,Vcmax25,Jmax25,GRNC,Einput,Eio);
%Result=Leaf(0.6,28,5,400,0,0,1)

Cii(i)=LeafResult(1);
NetAssimilationi(i)=LeafResult(2);
Gsi(i)=LeafResult(3);
LeafTemperaturei(i)=LeafResult(4);
Transpirationi(i)=LeafResult(5);
end
Air_CO2=0;
for j=1:40
WeatherTemp=25;
Air_CO2=10+40*j;
WeatherRH=0.6;
WeatherWind=5;
Radiation_PAR=470/2000*2000*0.85*0.85;%%%%
Radiation_NIR=0;
Radiation_LW=0;
PhotosynthesisType=1.1;
Vcmax25=100;
Jmax25=200;
GRNC=1;
LeafResult2=Leaf(WeatherRH,WeatherTemp,Air_CO2,WeatherWind,Radiation_PAR,Radiation_NIR,Radiation_LW,PhotosynthesisType,Vcmax25,Jmax25,GRNC,Einput,Eio);
%Result=Leaf(0.6,28,5,400,0,0,1)

HCii(j)=LeafResult2(1);
HNetAssimilationi(j)=LeafResult2(2);
HGsi(j)=LeafResult2(3);
HLeafTemperaturei(j)=LeafResult2(4);
HTranspirationi(j)=LeafResult2(5);
end



figure;
subplot(2,2,1); plot(Air_CO2i,NetAssimilationi,'k');hold on; plot(Air_CO2i,HNetAssimilationi,'r');xlabel('Air CO2');ylabel('Net CO2 assimilation (\mumol m^-^2 s^-^1)');
subplot(2,2,2); plot(Air_CO2i,Transpirationi,'k');hold on; plot(Air_CO2i,HTranspirationi,'r');xlabel('Air CO2');ylabel('Transpiration (\mumol m^-^2 s^-^1)');
subplot(2,2,3); plot(Air_CO2i,Gsi,'k');hold on; plot(Air_CO2i,HGsi,'r');xlabel('Air CO2');ylabel('Gs (mol m^-^2 s^-^1)');
subplot(2,2,4); plot(Air_CO2i,LeafTemperaturei,'k');hold on; plot(Air_CO2i,HLeafTemperaturei,'r');xlabel('Air CO2');ylabel('Leaf Temperature (^oC)');

R(:,1)=Air_CO2i';
R(:,2)=Cii';
R(:,3)=NetAssimilationi';
R(:,4)=Transpirationi';
R(:,5)=Gsi';
R(:,6)=LeafTemperaturei';
R(:,7)=HCii';
R(:,8)=HNetAssimilationi';
R(:,9)=HTranspirationi';
R(:,10)=HGsi';
R(:,11)=HLeafTemperaturei';


figure;
subplot(2,2,1); plot(R(:,2),R(:,3),'k');hold on; plot(R(:,7),R(:,8),'r');xlabel('Ci (\mumol mol^-^1)');ylabel('Net CO2 assimilation (\mumol m^-^2 s^-^1)');
subplot(2,2,2); plot(R(:,2),R(:,4),'k');hold on; plot(R(:,7),R(:,9),'r');xlabel('Ci (\mumol mol^-^1)');ylabel('Transpiration (\mumol m^-^2 s^-^1)');
subplot(2,2,3); plot(R(:,2),R(:,5),'k');hold on; plot(R(:,7),R(:,10),'r');xlabel('Ci (\mumol mol^-^1)');ylabel('Gs (mol m^-^2 s^-^1)');
subplot(2,2,4); plot(R(:,2),R(:,6),'k');hold on; plot(R(:,7),R(:,11),'r');xlabel('Ci (\mumol mol^-^1)');ylabel('Leaf Temperature (^oC)');


