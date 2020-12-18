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


function suc = PSTitle(m,p,n)
    if n==1
            if m ==1
                title('RUBP')
            elseif m==2
                title('PGA')
            elseif m==3
                title('DPGA')
            elseif m==4
                title('T3P')
            elseif m==5
                title('ADPG')
            elseif m==6
                title('FBP');
            elseif m==7
                title('E4P');
            elseif m==8
                title('S7P');
            elseif m==9
                title('SBP');
            elseif m==10
                title('ATP');
            elseif m==11
                title('NADPH');
            elseif m==12
                title('CO2');
            elseif m==13
                title('O2');
            elseif m==14
                title('HexP');
            elseif m==15
                title('PenP');
            else
                title('NONE');
            end
            suc = 1;
            xlabel('time(s)');
            ylabel('Concentration(mM)');
        elseif n ==2
            if m ==1
                title('Rubisco carboxylation')
            elseif m==2
                title('PGA kinase')
            elseif m==3
                title('GAPDH')
            elseif m==4
                title('T3P isomerase')
            elseif m==5
                title('Aldolase 1')
            elseif m==6
                title('FBPase');
            elseif m==7
                title('Transketolase 1');
            elseif m==8
                title('Aldolase 2');
            elseif m==9
                title('SBPase');
            elseif m==10
                title('Transketolase 2');
            elseif m==11
                title('Ru5P kinase');
            elseif m==12
                title('ATP synthase');
            elseif m==13
                title('ADPG Pyrophosphatase');
            elseif m==14
                title('DHAP export');
            elseif m==15
                title('PGA export');
            elseif m==16
                title('GAP export');
            else
                title('Starch synthase')
            end
                suc = 1;
                xlabel('time(s)');
                ylabel('Velocity(umol m-2 s-1)');
        end