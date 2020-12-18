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

RuACT_Con = RuACT_Ini(Begin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Getting the Velocity from RuACT_Rate.m %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global RuACT_OLD_TIME;
global RuACT_TIME_N;
global RuACT_VEL;
global RuACT_CON;

RuACT_OLD_TIME = 0;
RuACT_TIME_N = 1;

RuACT_VEL = zeros(1,3);    % Clean memory
RuACT_CON = zeros(3,1);    % Clean memory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign the major parameters in the system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

va1 = 1000;
n = 1;          % n can be sued to transfer certain variable. In this module is useless now. 

RuACT_Param = zeros(2,1);
RuACT_Param(1) = va1;
RuACT_Param(2) = n;

ModelComb = IniModelCom;        % Initialize the structure of the model, i.e. Is this model separate or combined with others. 

[Tt,d] = ode15s(@RuACT_MB,[0,time],RuACT_Con,options1,RuACT_Param);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get the output                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % Out put the Velocity Calculation

done = RuACT_Graph(Tt, d);
clock;

IniModelCom;