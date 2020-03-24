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

% %%%%%%%%%%%%%%%%%%%%%%%%
%   Getting the flux    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
global CO2A;
CO2A = zeros(5,2);
global PS_VEL;
global PS_OLD_TIME;
global PS_TIME_N;

iniModelCom;

global FIBF_PSPR_com; 
FIBF_PSPR_com = 0;


PS_OLD_TIME = 0;
PS_TIME_N = 1;
PS_VEL = zeros(17,1);   


Begin = 1;
PSs = PSInitial(Begin);

% %%%%%%%%%%%%%%%%%%%%%%%%%
%   CALCULATION            *
%%%%%%%%%%%%%%%%%%%%%%%%%%%

PS_Param = zeros(2,1);

ModelComb = IniModelCom;        % Initialize the structure of the model, i.e. Is this model separate or combined with others. 

PS_Param(1) = 1;            %

[Tt,d] = ode15s(@PSmb,[0,time],PSs,options1,PS_Param);

% %%%%%%%%%%%%%%%%%%%%%%%%%
%   OUTPUT                 *
%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUC = PS_OUT(Tt,d);

PSVEL = PS_VEL';
IniModelCom;