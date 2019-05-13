%   Copyright   Xin-Guang Zhu and Stephen P. Long, University of Illinois 
%   Copyright ? 2007

%   This file is part of CarbonMetabolism.

%    CarbonMetabolism is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.

%    CarbonMetabolism is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License (GPL)
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


% This is the function to output the different graphs from the simulation

function OUTSUC = CM_OUT(Tt, d)
global outA;
global outM;
global PS_VEL;
global PR_VEL;
global TandArateout;
global Arate;
% Get the concentration 
% figure;
% for m = 1:12
%     subplot(3,4,m);p = plot(Tt,d(:,m),'.');suc = PS_PR_Titl(m,p,1);
% end
% 
% figure;
% for m = 13:23
%     subplot(3,4,m-12);p = plot(Tt,d(:,m),'.');suc = PS_PR_Titl(m,p,1);
% end
% figure;
% for m = 24:35
%     subplot(3,4,m-23);p = plot(Tt,d(:,m),'.');ylabel('mM');xlabel('second'); suc = SUCS_AddTitle(m-23,p,1);
% end
TandCon(:,1)=Tt;
TandCon(:,2:36)=d;
[row,col]=size(TandCon);
% % TandConB=dataset({TandCon,'t','RuBP','PGA','DPGA','T3P','FBP','E4P','S7P','SBP','ATP','NADPH','CO2','O2','HexP','PenP','GCEA','GCA','PGCA','GCAc','GOAc','SERc','GLYc','HPRc','GCEAc','T3Pc','FBPc','HexPc','F26BPc','ATPc','ADPc','OPOPc','UDPGc','UTPc','SUCP','SUC','PGAc'});
% % BB=['t','RuBP','PGA','DPGA','T3P','FBP','E4P','S7P','SBP','ATP','NADPH','CO2','O2','HexP','PenP','GCEA','GCA','PGCA','GCAc','GOAc','SERc','GLYc','HPRc','GCEAc','T3Pc','FBPc','HexPc','F26BPc','ATPc','ADPc','OPOPc','UDPGc','UTPc','SUCP','SUC','PGAc'];
% % fid=fopen([BB,'C_Meta2.txt'],'w');
% % fprintf(fid,char(s));
% % fclose('all');
% % fid2=fopen([TandConB,'C_Meta2.txt'],'a');
% % fprintf(fid2, 'delimiter','\t','precision', '%.3f');
% % fclose('all');
% dlmwrite('C_Meta2.txt',TandConB,'delimiter','\t','precision', '%.3f');
% fid = fopen(outM,'wt'); 
% %fid = fopen('D6_C_Meta.txt','wt');  
% fprintf(fid,'%s','t RuBP PGA DPGA T3P FBP E4P S7P SBP ATP NADPH CO2 O2 HexP PenP GCEA GCA PGCA GCAc GOAc SERc GLYc HPRc GCEAc T3Pc FBPc HexPc F26BPc ATPc ADPc OPOPc UDPGc UTPc SUCP SUC PGAc');
% fprintf(fid,'%c\n',' ');    %»»ÐÐ
% for k=1:row; 
%     for m=1:36
%         p=num2str(TandCon(k,m));
%         fprintf(fid,'%s ',p);  
%     end
%     fprintf(fid,'%c\n',' ');  
% end
% fclose(fid); 


% Get the fluxes for the reactions in photorespiration pathway
CoeffVol = 30;
global VolRatioStCyto;      % VolRatioStCycto is the ratio between stroma and cytosol.
if VolRatioStCyto ==1
    for m=4:11
      PR_VEL(:,m) = PR_VEL(:,m)*CoeffVol;
    end
else
    PR_VEL(:,2) = PR_VEL(:,2)*CoeffVol;
    PR_VEL(:,3) = PR_VEL(:,3)*CoeffVol;
    PR_VEL(:,4) = PR_VEL(:,4)*CoeffVol;
    PR_VEL(:,5) = PR_VEL(:,5)*CoeffVol;
    PR_VEL(:,6) = PR_VEL(:,6)*CoeffVol;
    PR_VEL(:,7) = PR_VEL(:,7)*CoeffVol;
    PR_VEL(:,8) = PR_VEL(:,8)*CoeffVol;
    PR_VEL(:,9) = PR_VEL(:,9)*CoeffVol;
    PR_VEL(:,10) = PR_VEL(:,10)*CoeffVol;
    PR_VEL(:,11) = PR_VEL(:,11)*CoeffVol;
end

% figure;
%     for m = 1:10
%         subplot(3,4,m); p = plot(PR_VEL(:,1),PR_VEL(:,m+1),'.'); suc = PRespTitle(m,p,2);xlabel('time (s)');ylabel('uM m-2 s-1');
%     end

% Get the fluxes for reactions in photosynthesis pathway

PS_VEL = PS_VEL*CoeffVol; 
PS_VEL(1,:)=PS_VEL(1,:)/CoeffVol;

% figure;
% for m = 1:16
%     subplot(4,4,m); p = plot(PS_VEL(1,:), PS_VEL(1+m,:),'.'); suc = PSTitle(m,p,2);
% end

% The following plot the net carbon assimilation rate based on PS_VEL and PR_VEL.

CarbonRate = PS_VEL(2,:);
CO2Release = PR_VEL(:,9);
NetCarbon = CarbonRate' - CO2Release;
% figure;
% plot(PS_VEL(1,:),NetCarbon,'.');
% title('The net rate of carbon assimilation rate');
% xlabel('time(second)');
% ylabel('umol m-2 s-1');
TandArate(1,:)=PS_VEL(1,:);
TandArate(2,:)=NetCarbon;
TandArate=TandArate';
TandArateout=TandArate;
%dlmwrite('A_Meta2.txt',TandArate,'delimiter',' ','precision', '%.3f');
[row2,col2]=size(TandArate);
Arate=TandArate(row2,2);
% fid = fopen(outA,'wt');
% %fid = fopen('D7_A_Meta.txt','wt');  
% fprintf(fid,'%s','t PhotosynthesisRate');
% fprintf(fid,'%c\n',' ');   
% 
% for k=1:row2; 
%     for m=1:2
%         p=num2str(TandArate(k,m));
%         fprintf(fid,'%s ',p);   
%     end
%     fprintf(fid,'%c\n',' ');   
% end
% fclose(fid); 


% PSVEL = PS_VEL';
% PRVEL = PR_VEL';

global NetCarbon_plot;
NetCarbon_plot = NetCarbon;

global PS_VEL_gpmain;
PS_VEL_gpmain = PS_VEL;

PS_VEL(2:17,:) = PS_VEL(2:17,:)/CoeffVol;
PR_VEL(:,2:11) = PR_VEL(:,2:11)/CoeffVol;
OUTSUC = 1;
