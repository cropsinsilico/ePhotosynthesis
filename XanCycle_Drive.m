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


XanCycle_Con = XanCycle_Ini;


Begin = 1;
fin = SYSInitial(Begin);
global options1;
global tglobal;
time = tglobal;

Va1 = 1000;
n = 1;       

XanCycle_Param = zeros(2,1);
XanCycle_Param(1) = Va1;
XanCycle_Param(2) = n;

ModelComb = IniModelCom;         

[Tt,d] = ode15s(@XanCycle_MB,[0,time],XanCycle_Con,options1,XanCycle_Param);
 
    
for m = 1:4
    subplot(1,4,m);p = plot(Tt,d(:,m),'.');ylabel('mM');xlabel(' second');suc = XanCycle_AddTitle(m,p,1);
end
