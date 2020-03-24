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

function suc = FI_AddTitle(m,n,r)

    if r ==1
            if m ==1
                title('ISPr')
            elseif m==2
                title('cytf ox')
            elseif m==3
                title('ISP ox')
            elseif m==4
                title('ISPo-QH2')
            elseif m==5
                title('QH semi')
            elseif m==6
                title('cytbL ox')
           elseif m==7
                title('Q bound to Qi site')
            elseif m==8
                title('free Q in thylakoid memb')
            elseif m==9
                title('cytbH ox')
            elseif m==10
                title('Q-')
                
            elseif m==11
                title('Q2-')
            elseif m==12
                title('QH2')
            elseif m==13
                title('PC ox')
            elseif m==14
                title('P700')

            elseif m==16
                title('Pi')
            elseif m==17
                title('ATP')
            elseif m==18
                title('Ks')
            elseif m==19
                title('Mgs')
            elseif m==20
                title('Cls')
                
            elseif m==21
                title('Aip')
            elseif m==22
                title('U')
            elseif m==23
                title('A red')
            elseif m==24
                title('Fd red')
            elseif m==25
                title('BFHs + Hfs')
            elseif m==26
                title('BFHl + Hfl')
            elseif m==27
                title('PHs')
            elseif m==28
                title('PHl')
            elseif m==29
                title('NADPH')
            else
                title('None')
            end
            
        elseif r==2

            if m ==1
                title('Vbf1')
            elseif m==2
                title('Vbf2')
            elseif m==3
                title('Vbf3')
            elseif m==4
                title('Vbf4')
            elseif m==5
                title('Vbf5')
            elseif m==6
                title('Vbf6')
                
            elseif m==7
                title('Vbf7')
            elseif m==8
                title('Vbf8')
            elseif m==9
                title('Vbf9')
            elseif m==10
                title('Vbf10')
            elseif m==11
                title('Vbf11')
            elseif m==12
                title('Vqi')
                
            elseif m==13
                title('Vipc ')
            elseif m==14
                title('Vicp')
            elseif m==15
                title('Vinc')
            elseif m==16
                title('Vinp')
            elseif m==17
                title('Vdp')
                
            elseif m==18
                title('Vdc ')
            elseif m==19
                title('Vfp')
            elseif m==20
                title('Vfc')
            elseif m==21
                title('Vsfd')
            elseif m==22
                title('VsATP')
                
            elseif m==23
                title('VgPQH2 ')
            elseif m==24
                title('Vbf12')
            elseif m==25
                title('Vbf13')
            elseif m==26
                title('Vbf14')
            elseif m==27
                title('Vbf15')
              
            elseif m==28
                title('Vbf16')
            elseif m==29
                title('NADPH formation')
                
            else 
                title('none');
            end
        elseif r==3
            if m==1
                title('Membrane Potential');
            elseif m==2
                title('PHs');
            end
        end
            suc = 1;

