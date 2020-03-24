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



% trDynaPS_Drive.m
% This part include the function to begin the simulation.

% The time information is set in a global variable called tglobal in SYSInitial. 
time1 = clock;

Begin = 1;
fin = SYSInitial(Begin);
global options1;
global tglobal;
time = tglobal;

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

global RedoxReg_RA_com;     % This is the connection between Redox and RA.
RedoxReg_RA_com = 1;        % This means that the connection is there. 

global XanCycle_BF_com;
XanCycle_BF_com = 1;

global RROEA_EPS_com;
RROEA_EPS_com = 1;

global trDynaPS_SUCS_com;
trDynaPS_SUCS_com = 1;

global CO2A;
CO2A = zeros(5,1);

FPSReg_Con = FPSReg_Ini;

[FI_Param, BF_Param, PS_PR_Param, EPS_Param, RuACT_Param, XanCycle_Param, RROEA_Param, RedoxReg_Param] = ParamSet;

SUCS_Param = zeros(2,1);
SUCS_Param(1) = 0;
SUCS_Param(2) = 1;

[Tt,d] = ode15s(@FPSReg_mb,[0,time],FPSReg_Con,options1,BF_Param, FI_Param, PS_PR_Param, RuACT_Param, RedoxReg_Param, XanCycle_Param, RROEA_Param, SUCS_Param);

GraphSuc = FPSReg_Graph(Tt,d); 

    ATPActive = 0;
    BF_FI_com = 0;
    PR_PS_com = 0;
    FIBF_PSPR_com = 0;    
    RuACT_EPS_com = 0;
    RedoxReg_RA_com = 0;
    XanCycle_BF_com = 0;
    RROEA_EPS_com = 0;
    
    d = 0;
    
    global BF_VEL;
    global FI_VEL;
    global PS_VEL;
time2 =    clock;
RunTime = time2 - time1;
    
    IniModelCom;