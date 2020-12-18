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



function RedoxReg_Con = RedoxReg_Ini

BEGIN = 1;

global RedoxReg_OLD_TIME;
global RedoxReg_TIME_N;
global RedoxReg_VEL;
global RedoxReg_CON;

RedoxReg_OLD_TIME = 0;
RedoxReg_TIME_N = 1;
RedoxReg_VEL = zeros(1,3);     
RedoxReg_CON = zeros(3,1);     

 
RA_Con = RA_Ini;
RedoxReg_Con = zeros(5,1);

for m = 1:92
    RedoxReg_Con (m) = RA_Con(m);
end

Thion = 0.25;     % This is a wild guess
RedoxReg_Con(93) = Thion;             %


global RedoxReg_VMAX6;
global RedoxReg_VMAX9;
global RedoxReg_VMAX13;
global RedoxReg_VMAX16;

global V6;
global V9;
global V13;
global V16;

RedoxReg_VMAX6 = V6;
RedoxReg_VMAX9 = V9;
RedoxReg_VMAX13 = V13;
RedoxReg_VMAX16 = V16;

global RedoxReg_MP;

RedoxReg_MP(1,1) = 1000;
RedoxReg_MP(1,2) = -0.3;
RedoxReg_MP(1,3) = 0.5;

RedoxReg_MP(2,1) = 6;             % FBPase
RedoxReg_MP(2,2) = -0.305;
RedoxReg_MP(2,3) = 0.5;

RedoxReg_MP(3,1) = 9;             % SBPase
RedoxReg_MP(3,2) = -0.3;      
RedoxReg_MP(3,3) = 0.5;

RedoxReg_MP(4,1) = 13;            % PRK
RedoxReg_MP(4,2) = -0.295;
RedoxReg_MP(4,3) = 0.5;

RedoxReg_MP(5,1) = 16;            % ATPase
RedoxReg_MP(5,2) = -0.28;
RedoxReg_MP(5,3) = 0.5;


global Thio_Oxidation;       
global Fd_Thio_ET ;          

Thio_Oxidation = 0.1;     
Fd_Thio_ET  = 500;        
 

global ThioT;
ThioT = 0.5;                   

global BF_Pool;             
global BF2RedoxReg_Fdt;
BF2RedoxReg_Fdt = BF_Pool(6);   
