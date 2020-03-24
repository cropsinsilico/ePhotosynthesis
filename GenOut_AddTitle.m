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

function suc = GenOut_AddTitle(Response, m)
            if m ==1
                title('A')
            elseif m==2
                title('PGA')
            elseif m==3
                title('DPGA')
            elseif m==4
                title('T3P')
             elseif m==5
                title('FBP')
            elseif m==6
                title('E4P')
            elseif m==7
                title('S7P')
            elseif m == 8
                title('SBP')   
            elseif m ==9
                title('ATP')
            elseif m ==10
                title('NADPH')
            elseif m ==11
                title('HexP')
            elseif m==12
                title('PenP')
            elseif m == 13
                title('Pi')   
            elseif m ==14
                title('ADP')
            elseif m ==15
                title('RuBP')
            elseif m ==16
                title('Gcea')
            elseif m==17
                title('Gca')
            elseif m == 18
                title('Pga')   
            elseif m ==19
                title('Pgca')
            elseif m ==20
                title('Gcac')
            elseif m ==21
                title('Goac')
            elseif m==22
                title('SERc')
            elseif m == 23
                title('GLYc')   
            elseif m ==24
                title('Hprc')
            elseif m ==25
                title('Gceac')
            elseif m ==26
                title('None')
            elseif m==27
                title('Fdn')
            elseif m == 28
                title('PHs')   
            elseif m ==29
                title('PHl')
            elseif m ==30
                title('T3Pc')
            elseif m ==31
                title('FBPc')
            elseif m==32
                title('HexPc')
            elseif m == 33
                title('F26BPc')   
            elseif m ==34
                title('ATPc')
            elseif m ==35
                title('ADPc')
            elseif m ==36
                title('OPOPc')
            elseif m==37
                title('UDPGc')
            elseif m == 38
                title('UTPc')   
            elseif m ==39
                title('SUCP')
            elseif m ==40
                title('SUC')
            elseif m ==41
                title('PGAc')
             elseif m == 42
                title('V')   
            elseif m ==43
                title('A')
            elseif m ==44
                title('Z')
            elseif m ==45
                title('ABA')
            elseif m ==46
                title('TurgorPressure')
            elseif m ==47
                title('Gs')
            elseif m ==48
                title('Posm')                
               
            else 
                title('NONE')
            end
    
            
        switch Response
        case ('CO2')
            xlabel('[CO2] ppm');
        case ('O2')
            xlabel('[CO2] %');
        case ('light')
            xlabel('light uE');
        case ('enz')
            xlabel('Rubisco Act. (mM s-1)');
        otherwise
            xlabel('wrong');
        end
        ylabel('mM');    

    suc = 1;
    