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


function suc = PS_OUT(Tt,d)
global PS_VEL; 

for m = 1:6
    subplot(2,3,m);p = plot(Tt,d(:,m),'.');suc = PSTitle(m,p,1);
end
    pause;
for m = 7:12
    subplot(2,3,m-6);p = plot(Tt,d(:,m),'.');suc = PSTitle(m,p,1);
end
    pause;
    
for m = 10:15
    subplot(2,3,m-9);p = plot(Tt,d(:,m),'.');suc = PSTitle(m,p,1);
end
    pause;
    
global AVR; % The unit conversion between area and volume

PSVCOEFF = AVR;

PS_VEL = PS_VEL * PSVCOEFF; 

PS_VEL(1,:)=PS_VEL(1,:)/(PSVCOEFF);

for m = 1:8
    subplot(2,4,m); p = plot(PS_VEL(1,:), PS_VEL(1+m,:),'.'); suc = PSTitle(m,p,2);
end

pause;

for m = 9:17
    subplot(3,3,m-8); p = plot(PS_VEL(1,:), PS_VEL(1+m,:),'.'); suc = PSTitle(m,p,2);
end

pause; 

subplot(1,1,1); plot(PS_VEL(1,:),PS_VEL(2,:),'.'); title('CO2 assimilation rate');
