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

function RedoxReg_Vel = RedoxReg_Rate(t, RedoxReg_Con, RedoxReg_Param)

global RedoxReg_MP;

global RedoxReg_VMAX6;
global RedoxReg_VMAX9;
global RedoxReg_VMAX13;
global RedoxReg_VMAX16;

global BF2RedoxReg_Fdt;
global ThioT;

global Redox2PS_V6;
global Redox2PS_V9;
global Redox2PS_V13;
global Redox2PS_V16;

Fdn = RedoxReg_Con(24);
Fd = BF2RedoxReg_Fdt - Fdn;

Thion = RedoxReg_Con(93);
Thio = ThioT - Thion;

   RedoxReg_MP(1,3) = Thion/ThioT;
   
   TEMP = RedoxReg_MP(1,3);
   
   global RROEA_EPS_com;
    if RROEA_EPS_com == 1
        TEMP = 0.5;    
    end
    
global trDynaPS2RedReg_cal
    
if trDynaPS2RedReg_cal == 1      
   RedP = RedoxReg_MP(1,2) - 0.03 * log10(TEMP/(1-TEMP));   
  
   for index = 2:5
       
        RedPercent = RedoxReg_MP(index,3);
        MPE = RedoxReg_MP(index,2);
        pr = fsolve(@RedoxReg_FPercent, RedPercent,optimset('Display','off'),RedP,MPE);
        RedoxReg_MP(index,3) = pr;   
 
        switch RedoxReg_MP(index,1)
            
        case 6
            Redox2PS_V6 =  RedoxReg_VMAX6 * RedoxReg_MP(index,3);
            %Redox2PS_V6 = RedoxReg_VMAX6;
        case 9
            Redox2PS_V9 =  RedoxReg_VMAX9 * RedoxReg_MP(index,3);
            %Redox2PS_V9 = RedoxReg_VMAX9;
        case 13
            Redox2PS_V13 = RedoxReg_VMAX13 *  RedoxReg_MP(index,3);
            %Redox2PS_V13 = RedoxReg_VMAX13;
        case 16
            Redox2PS_V16 = RedoxReg_VMAX16 * RedoxReg_MP(index,3);
            %Redox2PS_V16 = RedoxReg_VMAX16;
        end
    end
end 

global Thio_Oxidation;    
global Fd_Thio_ET ;      

Vred = Fdn * Fd_Thio_ET * Thio/ThioT;
Vox = Thion * Thio_Oxidation;


RedoxReg_Vel = zeros(2,1);
RedoxReg_Vel(1) = Vred;
RedoxReg_Vel(2) = Vox;

 

global RedoxReg_OLD_TIME;
global RedoxReg_TIME_N;
global RedoxReg_VEL;
global RedoxReg_CON;

if (RedoxReg_TIME_N ==0)
    RedoxReg_TIME_N = 1;
end

if (t > RedoxReg_TIME_N)
    RedoxReg_TIME_N = RedoxReg_TIME_N + 1;
    RedoxReg_OLD_TIME = t;
end

RedoxReg_VEL(1,RedoxReg_TIME_N) = t;
RedoxReg_VEL(2,RedoxReg_TIME_N) = Vred;
RedoxReg_VEL(3,RedoxReg_TIME_N) = Vox;