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


function suc = RuACT_AddTitle(m,n,r)
        if r ==1
          
            if m ==1
                title('GAPDH')
            elseif m==2
                title('FBPase')
            elseif m==3
                title('SBPase')
            elseif m==4
                title('PRK')
             elseif m==5
                title('ATPase')
            elseif m==6
                title('ATPGPP')
            elseif m==7
                title('MDH')
            elseif m == 8
                title('Thio')   
            elseif m ==9
                title('Fd')
            else
                title('Rubisco Activase')
            end
            
        elseif r ==2
            if m ==1
                title('Thioo')
            end
        end
        
    suc = 1;