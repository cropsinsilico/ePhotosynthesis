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

function CO2AR = CM_Drive(pop, currentPop)

Begin = 1;
global tglobal;     % The total running time
tglobal = 10000;
global options1 
options1 = odeset('RelTol',1e-4);
time = tglobal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global variables used for obtaining flux and concentration data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global PS_PR_OLDTIME;
global PS_PR_TIME_N;
global PS_PR_VEL;

PS_PR_OLDTIME = 0;
PS_PR_TIME_N = 1;
PS_PR_VEL = zeros(27,1);      

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

global SUCS_OLD_TIME;
global SUCS_TIME_N;
global SUCS_VEL;
global SUCS_CON;

SUCS_OLD_TIME = 0;
SUCS_TIME_N = 1;
SUCS_VEL = zeros(1,3);    % Clean memory
SUCS_CON = zeros(3,1);    % Clean memory
%%%%%%%%%%%%%%%%%%%%%%%%
%   Initialation step %
%%%%%%%%%%%%%%%%%%%%%%%%

Begin = 1;
PS_PRs = zeros(23,1);
CMs = CM_Ini(Begin);

ModelComb = IModelCom;        % Initialize the structure of the model, i.e. Is this model standalone or combined with others. 

global PR_PS_com;    % This is a variable indicating whether the PR model is actually need to be combined with the Calvin cycle model or not. If 1 then means combined; 0 means not. 
PR_PS_com = 1;

global ATPActive;
ATPActive = 0;

global PSPR_SUCS_com;    % This is a variable indicating whether the PSPR model is actually need to be combined with SUCS or not. If 1 then means combined; 0 means not. 
PSPR_SUCS_com = 1;

%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculation  step %
%%%%%%%%%%%%%%%%%%%%%%%%

CM_Param = 0;

global CO2A;
CO2A = zeros(5,1);

SUCS_Param = zeros(2,1);
SUCS_Param(1) = 1;
SUCS_Param(2) = 1;

PS_PR_Param = 0;

suc = AssignVelocity(pop, currentPop);

[Tt,d] = ode15s(@CM_mb,[0,time],CMs, options1,PS_PR_Param, SUCS_Param);

global d_plot;
d_plot=d;

global Tt_plot;
Tt_plot = Tt;

%%%%%%%%%%%%%%%%%%%%%%%
%   output  step. 
%%%%%%%%%%%%%%%%%%%%%%%

PSPR_SUCS_com = 0;
IModelCom;

CO2AR = TargetFunVal