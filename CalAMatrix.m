%clear all;
    Edata=importdata('Einput7.txt');
    Eio=Edata.data(:,1);
    PPFDi=2000*0.85*0.85;
    %Tempi=25;
    global pcfactor;
    pcfactor=1.025;
    
for i=1:15 %Temp
    i
    Tempi=20+i;
    for j=1:120   %Ci
        j
        CO2i=j*10;
    Ao(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input5_0.txt',CO2i,PPFDi,Tempi,0,0,Eio);
    Ag(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input5_0.txt',CO2i,PPFDi,Tempi,1,0,Eio);

    end
    save AmatrixVc120_3P1.mat;
end

clear all;
 Edata=importdata('Einput7vc100.txt');
 Eio=Edata.data(:,1);
 PPFDi=2000*0.85*0.85;
 
 global pcfactor;
 pcfactor=1.025;
    
 
 for i=1:15 %Temp
    i
    Tempi=20+i;
    for j=1:120   %Ci
        j
        CO2i=j*10;
    Ao(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input5_0.txt',CO2i,PPFDi,Tempi,0,0,Eio);
    Ag(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input5_0.txt',CO2i,PPFDi,Tempi,1,0,Eio);

    end
    save AmatrixVc100_3P1.mat;
end
    Edata=importdata('Einput7.txt');
    Eio=Edata.data(:,1);
    PPFDi=2000*0.85*0.85;
    %Tempi=25;
    global pcfactor;
    pcfactor=1.056;
    
for i=1:15 %Temp
    i
    Tempi=20+i;
    for j=1:120   %Ci
        j
        CO2i=j*10;
    Ao(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input5_0.txt',CO2i,PPFDi,Tempi,0,0,Eio);
    Ag(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input5_0.txt',CO2i,PPFDi,Tempi,1,0,Eio);

    end
    save AmatrixVc120_3P2.mat;
end

clear all;
 Edata=importdata('Einput7vc100.txt');
 Eio=Edata.data(:,1);
 PPFDi=2000*0.85*0.85;
 
 global pcfactor;
 pcfactor=1.056;
    
 
 for i=1:15 %Temp
    i
    Tempi=20+i;
    for j=1:120   %Ci
        j
        CO2i=j*10;
    Ao(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input5_0.txt',CO2i,PPFDi,Tempi,0,0,Eio);
    Ag(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input5_0.txt',CO2i,PPFDi,Tempi,1,0,Eio);

    end
    save AmatrixVc100_3P2.mat;
end
% 
%  PPFDi=1200;
%     %Tempi=25;
%     
% for i=1:15 %Temp
%     i
%     Tempi=20+i;
%     for j=1:120   %Ci
%         j
%         CO2i=j*10;
%     Ao(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input4_0.txt',CO2i,PPFDi,Tempi,0,0,Eio);
%     Ag(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input4_0.txt',CO2i,PPFDi,Tempi,1,0,Eio);
% 
%     end
%     save AmatrixVc120.mat;
% end
% 
% 
%     Edata=importdata('Einput7.txt');
%     Eio=Edata.data(:,1);
%     PPFDi=2000;
%     %Tempi=25;
%     
% for i=1:15 %Temp
%     i
%     Tempi=20+i;
%     for j=1:120   %Ci
%         j
%         CO2i=j*10;
%     Ao(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input4_0.txt',CO2i,PPFDi,Tempi,0,0,Eio);
%     Ag(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input4_0.txt',CO2i,PPFDi,Tempi,1,0,Eio);
% 
%     end
%     save AmatrixPPFD2000Vc120.mat;
% end
% 
% 
%  PPFDi=1200;
%     %Tempi=25;
%     
% for i=1:15 %Temp
%     i
%     Tempi=20+i;
%     for j=1:120   %Ci
%         j
%         CO2i=j*10;
%     Ao(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input4_0.txt',CO2i,PPFDi,Tempi,0,0,Eio);
%     Ag(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input4_0.txt',CO2i,PPFDi,Tempi,1,0,Eio);
% 
%     end
%     save AmatrixPPFD2000Vc120.mat;
% end
% 
% 

% 
% 
%  PPFDi=1200;
%     %Tempi=25;
%     
% for i=1:15 %Temp
%     i
%     Tempi=20+i;
%     for j=1:120   %Ci
%         j
%         CO2i=j*10;
%     Ao(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input4_0.txt',CO2i,PPFDi,Tempi,0,0,Eio);
%     Ag(i,j,PPFDi)=EPS_Drive_GRNs('MeM_input4_0.txt',CO2i,PPFDi,Tempi,1,0,Eio);
% 
%     end
%     save AmatrixPPFD2000Vc100.mat;
% end