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
                title('A')
            elseif m==2
                title('U')
            elseif m==3
                title('P680Pheo')
            elseif m==4
                title('P680pPheon')
                
            elseif m==5
                title('P680pPheo')
            elseif m==6
                title('P680Pheon')
                
           elseif m==7
                title('Tyrosine')
                
            elseif m==8
                title('S1T')
            elseif m==9
                title('S2T')
            elseif m==10
                title('S3T')
            elseif m==11
                title('S0T')
            elseif m==12
                title('S1T+')
                
            elseif m==13
                title('S2T+')
            elseif m==14
                title('S3T+')
            elseif m==15
                title('S0T+')
            elseif m==16
                title('QAQB')
            elseif m==17
                title('QAnQB')
                
            elseif m==18
                title('QAQBn')
            elseif m==19
                title('QAnQBn')
            elseif m==20
                title('QAQB2n')
            elseif m==21
                title('QAnQB2n')
            elseif m==22
                title('PQn')
                
            else
                title('NONE')
            end
        elseif r==2

            if m ==1
                title('vA_d')
            elseif m==2
                title('vA_f')
            elseif m==3
                title('vAU')
            elseif m==4
                title('vUA')
            elseif m==5
                title('vU_f')
            elseif m==6
                title('vU_d')
                
            elseif m==7
                title('v1')
            elseif m==8
                title('vr1')
            elseif m==9
                title('vS1S2')
            elseif m==10
                title('vS2S3')
            elseif m==11
                title('vS3S0')
            elseif m==12
                title('vS0S1')
                
            elseif m==13
                title('vz_1 ')
            elseif m==14
                title('v1z_1')
            elseif m==15
                title('v2z_1')
            elseif m==16
                title('v3z_1')
            elseif m==17
                title('v0z_1')
                
            elseif m==18
                title('vz_2 ')
            elseif m==19
                title('v1z_2')
            elseif m==20
                title('v2z_2')
            elseif m==21
                title('v3z_2')
            elseif m==22
                title('v0z_2')
                
            elseif m==23
                title('v1z ')
            elseif m==24
                title('v2z')
            elseif m==25
                title('v3z')
            elseif m==26
                title('v0z')
            elseif m==27
                title('vAB1')
              
                
            elseif m==28
                title('vBA1')
            elseif m==29
                title('vAB2')
            elseif m==30
                title('vBA2')
            elseif m==31
                title('v3')
            elseif m==32
                title('vr3')
                
            elseif m==33
                title('v3n')
            elseif m==34
                title('v3rn')
            elseif m==35
                title('PQH2 oxi')
            elseif m==36
                title('IC ')
            elseif m==37
                title('IA')

            elseif m==38
                title('v2_1')
            elseif m==39
                title('v2_2')
            elseif m==40
                title('v2_00_1')
            elseif m==41
                title('v2_01_1 ')
            elseif m==42
                title('v2_02_1')
                
            elseif m==43
                title('v2_00_2')
            elseif m==44
                title('v2_01_2')
            elseif m==45
                title('v2_02_2')
            elseif m==46
                title('vr2_00_1 ')
            elseif m==47
                title('vr2_01_1')
                
            elseif m==48
                title('vr2_00_2')
            elseif m==49
                title('vr2_1')
            elseif m==50
                title('vr2_00_2')
            elseif m==51
                title('vr2_01_2 ')
            elseif m==52
                title('vr2_02_2')

            elseif m==52
                title('vr2_2')
            else 
                title('none');
                
            end
            
            elseif r==3
            if m==1
                title('fluo')
            elseif m==2
                title('Percent of QA red')
            elseif m==3
                title('heat dissipation')
            elseif m==4
                title('O2 evolution');
            elseif m==5
                title('PSII quantum yield')
           end
        
        end
        
            suc = 1;

