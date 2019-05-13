clear all;
% PCfactorMri=[0.964,1.002,0.960,0.984,0.986,1.011,0.972,0.999];
% inE=importdata('MeM_input11.txt');

% PCfactorMri=[0.700,0.687,1.01,1.037,0.961,0.980];
% inE=importdata('MeM_input9.txt');
PCfactorMri=[1,0.973,0.981,1.008];
inE=importdata('MeM_input12.txt');


[row,col]=size(inE.data);
% PCfactorMri=[0.973,1];

global pcfactor;  
% ProteinTotalRatio=0.973;
% pcfactor=1/ProteinTotalRatio;
% inE=importdata('MeM_input5_0.txt');
% Einput=inE.data;
Edata=importdata('Einput7.txt');
Eio=Edata.data(:,1);

for TFno=3:4
    TFno
ProteinTotalRatio=PCfactorMri(TFno);
pcfactor=1/ProteinTotalRatio;
Einput=inE.data(:,TFno);
for j=1:20
WeatherTemp=25;
Air_CO2=750;
WeatherRH=0.6;
WeatherWind=5;
Radiation_PAR=20*j-10;
Radiation_PARi(j)=Radiation_PAR;
Radiation_PARa=Radiation_PAR*0.85*0.85;
Radiation_NIR=0;
Radiation_LW=0;
PhotosynthesisType=1.1;
Vcmax25=100;
Jmax25=200;
GRNC=1;
LeafResult2=Leaf(WeatherRH,WeatherTemp,Air_CO2,WeatherWind,Radiation_PAR,Radiation_NIR,Radiation_LW,PhotosynthesisType,Vcmax25,Jmax25,GRNC,Einput,Eio);
%Result=Leaf(0.6,28,5,400,0,0,1)
HCii(j,TFno)=LeafResult2(1);
HNetAssimilationi(j,TFno)=LeafResult2(2);
HGsi(j,TFno)=LeafResult2(3);
HLeafTemperaturei(j,TFno)=LeafResult2(4);
HTranspirationi(j,TFno)=LeafResult2(5);
end
end

Radiation_PARi=Radiation_PARi*1E6/(2.35E5); %Convert W m^{-2} to u moles m^{-2} s^{-1}
% figure;
% subplot(2,2,1); plot(Radiation_PARi,NetAssimilationi,'r');hold on; plot(Radiation_PARi,HNetAssimilationi,'k');xlabel('PAR');ylabel('Net CO2 assimilation');
% subplot(2,2,2); plot(Radiation_PARi,Transpirationi,'r');hold on; plot(Radiation_PARi,HTranspirationi,'k');xlabel('PAR');ylabel('Transpiration');
% subplot(2,2,3); plot(Radiation_PARi,Gsi,'r');hold on; plot(Radiation_PARi,HGsi,'k'); xlabel('PAR');ylabel('Gs');
% subplot(2,2,4); plot(Radiation_PARi,LeafTemperaturei,'r');hold on; plot(Radiation_PARi,HLeafTemperaturei,'k');xlabel('PAR');ylabel('Leaf Temperature');

figure;
plot(Radiation_PARi,HNetAssimilationi(:,1),'k');hold on;
plot(Radiation_PARi,HNetAssimilationi(:,6),'r');
plot(Radiation_PARi,HNetAssimilationi(:,8),'b');
plot(Radiation_PARi,HNetAssimilationi(:,1),'k.');
plot(Radiation_PARi,HNetAssimilationi(:,5),'r.');
plot(Radiation_PARi,HNetAssimilationi(:,7),'b.');
xlabel('PAR (\mumol mol^-^1)');
ylabel('CO_2 uptale rate (\mu mol m^-^2 s^-^1)');
xlim([0,1500]);
ylim([0,30]);


figure;
subplot(1,3,1);plot(AA(:,1),AA(:,2),'k');hold on;plot(AA(:,1),AA(:,3),'k.');
subplot(1,3,2);plot(AA(:,1),AA(:,4),'k');hold on;plot(AA(:,1),AA(:,5),'k.');
subplot(1,3,3);plot(AA(:,1),AA(:,6),'k');hold on;plot(AA(:,1),AA(:,7),'k.');

plot(AA(:,1),HNetAssimilationi(:,6),'r');
plot(AA(:,1),HNetAssimilationi(:,8),'b');
plot(AA(:,1),HNetAssimilationi(:,1),'k.');
plot(AA(:,1),HNetAssimilationi(:,5),'r.');
plot(AA(:,1),HNetAssimilationi(:,7),'b.');
xlabel('PAR (\mumol mol^-^1)');
ylabel('CO_2 uptale rate (\mu mol m^-^2 s^-^1)');
xlim([0,1500]);
ylim([0,30]);



