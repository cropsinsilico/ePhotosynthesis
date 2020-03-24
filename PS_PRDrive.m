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



clear ; 


Begin = 1;
fin = SYSInitial(Begin);
global options1;
global tglobal;
time = tglobal;


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

%%%%%%%%%%%%%%%%%%%%%%%%
%   Initialation step %
%%%%%%%%%%%%%%%%%%%%%%%%
Begin = 1;

PS_PRs = zeros(24,1);

ModelComb = IniModelCom;       

global PR_PS_com;     
PR_PS_com = 1;

PS_PRs = PS_PRIni(Begin);

PS_PR_Param = 0;

global CO2A;
CO2A = zeros(5,1);

[Tt,d] = ode15s(@PS_PRmb,[0,time],PS_PRs,options1,PS_PR_Param);


success = PSPR_Out(Tt, d);


PR_PS_com = 0;
IniModelCom;