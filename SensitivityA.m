clear all;
        k=2000*0.85*0.85;
        PPFDi=k;
        Tempi=25;
        CO2iH=550*0.7;
        CO2iL=370*0.7;
        Edata=importdata('Einput7.txt');
        Eio=Edata.data(:,1);
        
        global pcfactor;
         pcfactor=1.045;%1/0.973;%1.423;
         Eiot=Eio;
% %          Eiot(1:27)=Eio(1:27)*pcfactor;
% %          Eiot(38)=Eio(38)*pcfactor;
% %          Eiot(45)=Eio(45)*pcfactor;
% %      if pcfactor==1.059;
% %          Eiot(1)=Eiot(1)/pcfactor;
% %      end
     if pcfactor==1.045;%1.423;
         Eiot(1)=Eio(1)*pcfactor;
         Eiot(4)=Eio(4)*pcfactor;
         Eiot(5)=Eio(5)*pcfactor;
         Eiot(8)=Eio(8)*pcfactor;
         Eiot(11)=Eio(11)*pcfactor;
         Eiot(21)=Eio(21)*pcfactor;
         Eiot(22)=Eio(22)*pcfactor;
         Eiot(25)=Eio(25)*pcfactor;
         Eiot(38)=Eio(38)*pcfactor;
     end
      pcfactor=1;
     
%         if pcfactor==1.423;
%          Eiot(4)=Eio(4)*pcfactor;
%          Eiot(5)=Eio(5)*pcfactor;
%          Eiot(8)=Eio(8)*pcfactor;
%          Eiot(11)=Eio(11)*pcfactor;
%          Eiot(21)=Eio(21)*pcfactor;
%          Eiot(22)=Eio(22)*pcfactor;
%          Eiot(25)=Eio(25)*pcfactor;
%          Eiot(38)=Eio(38)*pcfactor;
%      end

inE=importdata('MeM_input5_0.txt');
Einput=inE.data;
        
        StanAo=EPS_Drive_GRNs(Einput,CO2iL,PPFDi,Tempi,0,0,Eio);
        StanAe=EPS_Drive_GRNs(Einput,CO2iH,PPFDi,Tempi,0,0,Eio);
        StanAg=EPS_Drive_GRNs(Einput,CO2iH,PPFDi,Tempi,1,0,Eiot);
        
        
        for j=1:66   %Ci
        j
        ENum=j;
        Eic=Eio;
        Eict=Eiot;
        Eic(j)=Eio(j)*0.99;
        Eict(j)=Eiot(j)*0.99;
        %Ao(j)=EPS_Drive_GRNs('MeM_input5_0.txt',CO2iL,PPFDi,Tempi,0,0,Eic);
        Ae(j)=EPS_Drive_GRNs(Einput,CO2iH,PPFDi,Tempi,0,0,Eic);
        Ag(j)=EPS_Drive_GRNs(Einput,CO2iH,PPFDi,Tempi,1,0,Eict);
        Eid=Eio;
        Eidt=Eiot;
        Eid(j)=Eio(j)*1.01;
        Eidt(j)=Eiot(j)*1.01;
        %Aod(j)=EPS_Drive_GRNs('MeM_input5_0.txt',CO2iL,PPFDi,Tempi,0,0,Eid);
        Aed(j)=EPS_Drive_GRNs(Einput,CO2iH,PPFDi,Tempi,0,0,Eid);
        Agd(j)=EPS_Drive_GRNs(Einput,CO2iH,PPFDi,Tempi,1,0,Eidt);        
        
        end
        %Seno=(Aod-Ao)/StanAo/0.02;
        Sene=(Aed-Ae)/StanAe/0.02;
        Seng=(Agd-Ag)/StanAg/0.02;
        
%         Seno=(StanAo-Ao)/StanAo/0.01;
%         Sene=(StanAe-Ae)/StanAe/0.01;
%         Seng=(StanAg-Ag)/StanAe/0.01;
        T(1:66,1)=Sene';
        T(1:66,2)=Sene';
        T(1:66,3)=Seng';
        T(67,1)=StanAo;
        T(67,2)=StanAe;
        T(67,3)=StanAg;