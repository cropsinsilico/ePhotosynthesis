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

function EPS_Con = EPS_Ini(BEGIN)

Begin = 1;

FIBFsuc = FIBF_Ini(Begin);
FI_Con = FI_Ini(Begin);
BF_Con = BF_Ini(Begin);

%% Step 1 Get the initialization of the variables for BF

global BF_OLD_TIME;
global BF_TIME_N;
global BF_VEL;
global BF_CON;

BF_OLD_TIME = 0;
BF_TIME_N = 1;

BF_VEL = zeros(1,5);    % Clean memory
BF_CON = zeros(1,5);    % Clean memory

%% Get the initialization of the variables for FI

global FI_OLD_TIME;
global FI_TIME_N;
global FI_VEL;
global FI_CON;

FI_OLD_TIME = 0;
FI_TIME_N = 1;

FI_VEL = zeros(1,5);    % Clean memory
FI_CON = zeros(5,1);    % Clean memory


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
    
FIBF_Con(52) = 10^8 * 0.5;    % The initialization of the initial rate constant for heat dissipation
    
% Second, try to get the iniitalzation files for the PSPR model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global variables used for obtaining flux and concentration data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global PS_PR_OLDTIME;
global PS_PR_TIME_N;
global PS_PR_VEL;

PS_PR_OLDTIME = 0;
PS_PR_TIME_N = 1;
PS_PR_VEL = zeros(27,1);        % Store the flux value

global PS_OLD_TIME;
global PS_TIME_N;
global PS_VEL;
PS_OLD_TIME = 0;
PS_TIME_N= 0;
PS_VEL = zeros(1,1);

global PR_OLD_TIME;
global PR_TIME_N;
global PR_VEL;
PR_OLD_TIME = 0;
PR_TIME_N = 1;
PR_VEL = zeros(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%
%   Initialation step %
%%%%%%%%%%%%%%%%%%%%%%%%

Begin = 1;
CMs = CM_Ini(Begin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Use the initialized variables to construct variables that will be transfered to the Drive file. %%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for m = 1:52
    EPS_Con(m) = FIBF_Con(m);
end

for m = 1:35
    EPS_Con(m+52) = CMs(m);
end

