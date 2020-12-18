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

% EPS_Drive.m
% This part include the function to begin the simulation.

% The time information is set in a global variable called tglobal in SYSInitial. 
clear; 

Begin = 1;
fin = SYSInitial(Begin);
global options1;
global tglobal;
time = tglobal;

%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculation  step  %
%%%%%%%%%%%%%%%%%%%%%%%%

global ATPActive;
ATPActive = 0;

global EPS_ATP_Rate;        % Indicate in the beginning there is no ATP synthesis activity.
EPS_ATP_Rate = 0;

ModelComb = IniModelCom;        % Initialize the structure of the model, i.e. Is this model separate or combined with others. 

global BF_FI_com;            % The combination of BF and FI model 
BF_FI_com = 1;

global PR_PS_com;    % This is a variable indicating whether the PR model is actually need to be combined with PS or not. If 1 then means combined; 0 means not. 
PR_PS_com = 1;

global FIBF_PSPR_com; % 1 means that the overall EPS model is used. 0 means partial model of FIBF is used. 
FIBF_PSPR_com = 1;    

global RuACT_EPS_com;     % A global variable to indicate whether the RuACT is run by itself or combined with others. 
RuACT_EPS_com = 1;        % Since this is run within this program, it is combinbed, therefore, it is assigned value 1, otherwise, assign value 0. 

global EPS_SUCS_com;
EPS_SUCS_com = 1;

global PSPR_SUCS_com;    % This is a variable indicating whether the PSPR model is actually need to be combined with SUCS or not. If 1 then means combined; 0 means not. 
PSPR_SUCS_com = 1;

global CO2A;
CO2A = zeros(5,1);

RA_Con = RA_Ini;

va1 = 0;
global PS12ratio; % The ratio of the PSI unit to the PSII unit
BF_Param = zeros(5,1);
BF_Param(1) = va1;
BF_Param(2) = PS12ratio;

FI_Param = zeros(5,1);
FI_Param(1) = va1;
FI_Param(2) = PS12ratio;

PS_PR_Param = 0;
EPS_Param = 0;

RuACT_Param = zeros(2,1);
RuACT_Param(1) = va1;
RuACT_Param(2) = PS12ratio;
Sucs_Param = 0; 


[Tt,d] = ode15s(@RA_mb,[0,time],RA_Con,options1,BF_Param, FI_Param, PS_PR_Param, RuACT_Param, Sucs_Param);

%   The following section deals with the data output of the program.

    M = size(d);
    a = M(1);
    
    FIBF = zeros(a,52);
    PSPR = zeros(a,23);
    
    for m = 1:52
         FIBF(:,m)= d(:,m);
    end
    
    for m = 1:24
         PSPR(:,m)=d(:,m+52);
    end
    
    global FI;
    global BF;
    
    global PS;
    global PR;
    
    OutPutLight = FIBF_Out(Tt, FIBF);
    OutputPSPR = PSPR_Out(Tt, PSPR);
    
    % output the rate of CO2 UTPAKE.
    global RuACT_VEL;
    global PR_VEL;
    CarbonRate = RuACT_VEL(:,6) * 37.5;
    OxygenationRate = RuACT_VEL(:,7) * 37.5;
    
    CO2Release = PR_VEL(:,9)*37.5 ;
    Assim = CarbonRate - CO2Release; 
    
    newplot;
    subplot(2,2,1)
    p = plot((PR_VEL(:,1)), Assim,'.');ylabel('Assimilation Rate');xlabel('second');title('Net Assimilation Rate');
    subplot(2,2,2)    
    p = plot((PR_VEL(:,1)), CarbonRate,'.');ylabel('Carboxylation rate');xlabel('second');title('RuBP Carboxylation rate');
    subplot(2,2,3)    
    p = plot((PR_VEL(:,1)), OxygenationRate,'.');ylabel('Oxygenation Rate');xlabel('second');title('Oxygenation Rate');
    subplot(2,2,4)    
    p = plot((PR_VEL(:,1)), CO2Release,'.');ylabel('CO2 release from GDC ');xlabel('second');title('CO2 release from GDC');
    pause;
    
    d_RuACT= d(:,89:92);
    done = RuACT_Graph(Tt, d_RuACT);
    
    ATPActive = 0;
    
    BF_FI_com = 0;
    PR_PS_com = 0;
    FIBF_PSPR_com = 0;    
    RuACT_EPS_com = 0;
    global BF_VEL;
    global FI_VEL;
    global PS_VEL;
    
    IniModelCom;