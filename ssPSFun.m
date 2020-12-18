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


function [ATP, w] =  ssPSFun (VcmaxT, JmaxT, temp, CO2, Light)

global kmCO2; 
global kmO2;
global O2;
global GammaStar;
global Rd;

done = ssPSIni(temp);

Ci = CO2;
wc = VcmaxT * (Ci - GammaStar)/(Ci + kmCO2 * (1 + O2/kmO2));
wj = JmaxT * (Ci - GammaStar)/(4.5 * Ci + 10.5 * GammaStar);
w = min(wc, wj);

Vm = 88.6 *10^(-3);
at = 12.6 * Vm;
p = 2.5 * Vm;
Vr = Vm * 2.27;

tC = at - p * wc/wj;
tJ = (at - p) * (Vr /Vm - 1)/(wc * Vr/(wj * Vm) - 1);

To = 1;      

if wc < wj
    ATP = To + tC;
else
    ATP = tJ;
end
                        
ATP = ATP/5 * 1.5;       
