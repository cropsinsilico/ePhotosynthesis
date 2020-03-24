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


Begin = 1;
fin = SYSInitial(Begin);
global options1;
global tglobal;
time = tglobal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get the initial condition %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FI_Con = FI_Ini(Begin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Getting the Velocity from FI_Rate.m %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global FI_OLD_TIME;
global FI_TIME_N;
global FI_VEL;
global FI_CON;

FI_OLD_TIME = 0;
FI_TIME_N = 1;

FI_VEL = zeros(1,5);    % Clean memory
FI_CON = zeros(5,1);    % Clean memory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign the major parameters in the system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


va1 = 0;
global PS12ratio;

FI_Param = zeros(5,1);
FI_Param(1) = va1;
FI_Param(2) = PS12ratio;

ModelComb = IniModelCom;        % Initialize the structure of the model, i.e. Is this model separate or combined with others. 

[Tt,d] = ode15s(@FI_MB,[0,time],FI_Con,options1,FI_Param);
    
    for m = 1:10
        subplot(2,5,m);p = plot(log10(Tt),d(:,m),'.');ylabel('micro mol/m2');xlabel('log 10 second');suc = FI_AddTitle(m,p,1);
    end
    pause;

    for m = 11:20
        subplot(2,5,m-10);p = plot(log10(Tt),d(:,m),'.');ylabel('micro mol/m2');xlabel('log 10 second');suc = FI_AddTitle(m,p,1);
    end
    pause;

    for m = 1:10
        subplot(2,5,m);p = plot(log10(FI_VEL(:,1)),FI_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,2);
    end
    pause;

    for m = 11:20
        subplot(2,5,m-10);p = plot(log10(FI_VEL(:,1)),FI_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,2);
    end
    pause;

    for m = 21:30
        subplot(2,5,m-20);p = plot(log10(FI_VEL(:,1)),FI_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,2);
    end
    pause;

    for m = 31:40
        subplot(2,5,m-30);p = plot(log10(FI_VEL(:,1)),FI_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,2);
    end
    pause;
    

    for m = 41:50
        subplot(2,5,m-40);p = plot(log10(FI_VEL(:,1)),FI_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,2);
    end
    pause;
    
    for m = 51:53
        subplot(1,3,m-50);p = plot(log10(FI_VEL(:,1)),FI_VEL(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,2);
    end
    pause;
    
    for m = 1:5
        subplot(1,5,m);p = plot(log10(FI_VEL(:,1)),FI_CON(:,m+1),'.');ylabel('micro mol/m2/s');xlabel('log 10 second');suc = FI_AddTitle(m,p,3);
    end
IniModelCom;