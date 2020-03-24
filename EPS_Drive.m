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



clear; 

Begin = 1;
fin = SYSInitial(Begin);
global options1;
global tglobal;
time = tglobal;

%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculation  step %
%%%%%%%%%%%%%%%%%%%%%%%%
EPS_Con = EPS_Ini;

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
Sucs_Param = 0; 
% Get the initial variables for SUCS.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       END                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a variable to indicate whether BF is run by itself or together. 

ModelComb = IniModelCom;        % Initialize the structure of the model, i.e. Is this model separate or combined with others. 

global BF_FI_com;
BF_FI_com = 1;

global PR_PS_com;    % This is a variable indicating whether the PR model is actually need to be combined with PS or not. If 1 then means combined; 0 means not. 
PR_PS_com = 1;

global FIBF_PSPR_com; % 1 means that the overall EPS model is used. 0 means partial model is used. 
FIBF_PSPR_com = 1;    

global ATPActive;
ATPActive = 0;

global EPS_ATP_Rate;        % Indicate in the beginning there is no ATP synthesis activity.
EPS_ATP_Rate = 0;

global EPS_SUCS_com;
EPS_SUCS_com = 1;

global PSPR_SUCS_com;    % This is a variable indicating whether the PSPR model is actually need to be combined with SUCS or not. If 1 then means combined; 0 means not. 
PSPR_SUCS_com = 1;

global CO2A;
CO2A = zeros(5,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Tt,d] = ode15s(@EPS_mb,[0,time],EPS_Con,options1,BF_Param, FI_Param, PS_PR_Param, Sucs_Param);

%   The following section deals with the data output of the program.
GraphSuc = EPS_Graph(Tt,d);    
% Here is the place to recover the status of the regulatory variables. 
    
%     ATPActive = 0;
%     
%     FIBF_PSPR_com = 0;
%     BF_FI_com = 0;
%     PR_PS_com = 0;
%     IniModelCom;
    
    global PS_VEL; 
    global PR_VEL;
    global PS_CON; 
    global FI_VEL; 
    global FI_CON; 
    global BF_VEL; 
    global BF_CON; 