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

FIBFsuc = FIBF_Ini(Begin);
FI_Con = FI_Ini(Begin);
BF_Con = BF_Ini(Begin);


global FIBF_Pool;
global FI_Pool;
global BF_Pool;
FI_Pool(2) = FIBF_Pool(1);
BF_Pool(8) = FIBF_Pool(1);

% Initial concentration for FIBF_Con
    for m = 1:29
        FIBF_Con(m) = BF_Con(m);
    end
    
    for m = 1:22
        FIBF_Con(m+29) = FI_Con(m);
    end
  
FIBF_Con(52) = 0.5 * 10^8;    % The initialization of the initial rate constant for heat dissipation
    
%--------------------------------------------------------------|
%--------------------------------------------------------------|
% Calculation step                                             |
%---------------------- ---------------------------------------|

%% Step 1 Get the initialization of the variables for BF

global BF_OLD_TIME;
global BF_TIME_N;
global BF_VEL;
global BF_CON;

BF_OLD_TIME = 0;
BF_TIME_N = 1;

BF_VEL = zeros(1,5);    % Clean memory
BF_CON = zeros(1,5);    % Clean memory

va1 = 0;
global PS12ratio; % The ratio of the PSI unit to the PSII unit
BF_Param = zeros(5,1);
BF_Param(1) = va1;
BF_Param(2) = PS12ratio;

FI_Param = zeros(5,1);
FI_Param(1) = va1;
FI_Param(2) = PS12ratio;

% Step 2 Get teh initialization of the variables for FI.
global FI_OLD_TIME;
global FI_TIME_N;
global FI_VEL;
global FI_CON;

FI_OLD_TIME = 0;
FI_TIME_N = 1;

FI_VEL = zeros(1,5);    % Clean memory
FI_CON = zeros(5,1);    % Clean memory

%%  Step 3 Define the variables for the FIBF

ModelComb = IniModelCom;        % Initialize the structure of the model, i.e. Is this model separate or combined with others. 

global BF_FI_com;            % The combination of BF and FI model 
BF_FI_com = 1;


[Tt,d] = ode15s(@FIBF_MB,[0,time],FIBF_Con,options1,BF_Param, FI_Param);                                         
    Success = FIBF_Out(Tt,d);
    
    
% Some of the parameters need to return to its original value. 
    BF_FI_com = 0;
IniModelCom;