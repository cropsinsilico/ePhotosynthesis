%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%   Copyright   Xin-Guang Zhu, Yu Wang, Donald R. ORT and Stephen P. LONG
%CAS-MPG Partner Institute for Computational Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
%China Institute of Genomic Biology and Department of Plant Biology, Shanghai Institutes for Biological Sciences, CAS, Shanghai,200031 
%University of Illinois at Urbana Champaign
%Global Change and Photosynthesis Research Unit, USDA/ARS, 1406 Institute of Genomic Biology, Urbana, IL 61801, USA.
 
%   This file is part of e-photosynthesis.
 
%    e-photosynthesis is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; 
 
%    e-photosynthesis is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
 
%    You should have received a copy of the GNU General Public License (GPL)
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% FI_Drive.m
Begin = 1;
fin = SYSInitial(Begin);
global options1;
global tglobal;
time = tglobal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get the initial condition %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BF_Con = BF_Ini(Begin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Getting the Velocity from BF_Rate.m %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global BF_OLD_TIME;
global BF_TIME_N;
global BF_VEL;
global BF_CON;

BF_OLD_TIME = 0;
BF_TIME_N = 1;

BF_VEL = zeros(1,5);    % Clean memory
BF_CON = zeros(1,5);    % Clean memory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign the major parameters in the system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

va1 = 0;  % It is defined now in the conditition routine. 
global PS12ratio; % The ratio of the PSI unit to the PSII unit

BF_Param = zeros(5,1);
BF_Param(1) = va1;
BF_Param(2) = PS12ratio;


ModelComb = IniModelCom;       

[Tt,d] = ode15s(@BF_MB,[0,time],BF_Con,options1,BF_Param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get the output                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Out put the concentrtaion

    for m = 1:10
        subplot(2,5,m);p = plot((Tt),d(:,m),'.');ylabel('micro mol/m2');xlabel('second');suc = BF_AddTitle(m,p,1);
    end
    pause;

    for m = 11:14
        subplot(1,4,m-10);p = plot((Tt),d(:,m),'.');ylabel('micro mol/m2');xlabel(' second');suc = BF_AddTitle(m,p,1);
    end
    pause;
   
    for m = 16:20
        subplot(1,5,m-15);p = plot((Tt),d(:,m),'.');ylabel('mmol/l');xlabel(' second');suc = BF_AddTitle(m,p,1);
    end
    pause;

    for m = 21:24
        subplot(2,2,m-20);p = plot((Tt),d(:,m),'.');ylabel('micromol/m2');xlabel('second');suc = BF_AddTitle(m,p,1);
    end
    pause;
    
    for m = 25:26
        subplot(1,2,m-24);p = plot((Tt),d(:,m),'.');ylabel('milimole/l');xlabel('second');suc = BF_AddTitle(m,p,1);
    end
    pause;
    
    for m = 27:28
        subplot(1,2,m-26);p = plot((Tt),d(:,m),'.');ylabel('PH');xlabel('second');suc = BF_AddTitle(m,p,1);
    end
    pause;
    
    newplot;
    p = plot((Tt),d(:,29),'.');ylabel('mmol l-1');xlabel('second');suc = BF_AddTitle(29,p,1);
    pause;
    
    % Output the Velocity
    for m = 1:10
        subplot(2,5,m);p = plot((BF_VEL(:,1)),BF_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('second');suc = BF_AddTitle(m,p,2);
    end
    pause;
    
    newplot;
    p = plot((BF_VEL(:,1)),BF_VEL(:,12),'.');ylabel('mmol/l/s');xlabel('second');suc = BF_AddTitle(11,p,2);
    pause;
    
    for m = 12:20
        subplot(2,5,m-11);p = plot((BF_VEL(:,1)),BF_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('second');suc = BF_AddTitle(m,p,2);
    end
    pause;

    for m = 21:22
        subplot(2,4,m-20);p = plot((BF_VEL(:,1)),BF_VEL(:,m+1),'.');ylabel('mmol/l/s');xlabel('second');suc = BF_AddTitle(m,p,2);
    end
    pause;

    newplot;
    p = plot((BF_VEL(:,1)),BF_VEL(:,23),'.');ylabel('mircomol/m2/s');xlabel('second');suc = BF_AddTitle(23,p,2);
    pause;
   
    for m = 24:26
        subplot(1,3,m-23);p = plot((BF_VEL(:,1)),BF_VEL(:,m+1),'.');ylabel('mmol/l/s');xlabel('second');suc = BF_AddTitle(m,p,2);
    end
    pause;
    
    for m = 27:28
        subplot(1,2,m-26);p = plot((BF_VEL(:,1)),BF_VEL(:,m+1),'.');ylabel('micromol /m2 leaf/s');xlabel('second');suc = BF_AddTitle(m,p,2);
    end
    pause;

    newplot;
    p = plot((BF_VEL(:,1)),BF_VEL(:,30),'.');ylabel('mmol/l/s');xlabel('second');suc = BF_AddTitle(29,p,2);
    pause;
    
    newplot;
    p = plot((BF_CON(:,1)),BF_CON(:,2),'.');ylabel('V(from stroma to lumen)');xlabel('second');suc = BF_AddTitle(1,p,3);
    
    IniModelCom;