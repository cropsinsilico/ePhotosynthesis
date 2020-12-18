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


function CA_Table(MetLength,VmaxLength,VmaxVary,VmaxLabels,MetConCoeff,MetaboliteLabels,OriginalValues,ControlCoefficientsTable,DriveSelector)

Table = cell(MetLength+4,VmaxLength+1);
for j=1:MetLength+4
    for h=1:VmaxLength+1
        Table{j,h} = ' ';
    end
end

count1 = 1;
count2 = 1; 
count3 = 1;
count4 = 1; 

for a=1:VmaxLength
    if VmaxVary(a,1) == 1
        Table{1,count1+1}  = VmaxLabels(a,:);
        count1 = count1 + 1;
    end
end

for s=1:MetLength
    if MetConCoeff(s,1) == 1
        Table{count2+2,1} = MetaboliteLabels(s,:);
        count2 = count2 + 1;
    end
end
Table{2,1}  = 'Original Values';

for n=1:VmaxLength
    if VmaxVary(n,1) == 1
        Table{2,count3+1} = OriginalValues(n,1); 
        count3 = count3 + 1;
    end
end


for g=1:MetLength
    if MetConCoeff(g,1) == 1
        count5 = 1;
        for d=1:VmaxLength
            if VmaxVary(d,1) == 1                             
                Table{count4+2,count5+1} = ControlCoefficientsTable(g,d);
                count5 = count5 + 1;        
            end
        end
        count4 = count4 + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Table{count2+4,1} = 'Drive Used';

if DriveSelector == 1
    Table{count2+4,2} = 'PS Drive';
elseif DriveSelector == 2
    Table{count2+4,2} = 'PS-PR Drive';
elseif DriveSelector == 3
    Table{count2+4,2} = 'trDynaPS Drive';
elseif DriveSelector == 4
    Table{count2+4,2} = 'CM_Drive';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlswrite('ControlCoefficients.xls',Table,'A1');
