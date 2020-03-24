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


% PSPR_Out.m    This program provide the output for the model of the PS and PR.

function Success = PSPR_Out(Tt, d)

global PS_VEL;
global PR_VEL;


% Get the concentration 
for m = 1:6
    subplot(2,3,m);p = plot(Tt,d(:,m),'.');suc = PS_PR_Titl(m,p,1);
end
pause;

for m = 7:12
    subplot(2,3,m-6);p = plot(Tt,d(:,m),'.');suc = PS_PR_Titl(m,p,1);
end
pause;

for m = 13:18
    subplot(2,3,m-12);p = plot(Tt,d(:,m),'.');suc = PS_PR_Titl(m,p,1);
end
pause;


for m = 19:24
    subplot(2,3,m-18);p = plot(Tt,d(:,m),'.');suc = PS_PR_Titl(m,p,1);
end
pause;

DOTEST = 0; 

if DOTEST ==0
% Get the fluxes for the reactions in photorespiration pathway      
	global AVR; 
	Coeff = AVR ;
	PR_VEL(:,2) = PR_VEL(:,2)*Coeff;
	PR_VEL(:,3) = PR_VEL(:,3)*Coeff;
	PR_VEL(:,4) = PR_VEL(:,4)*Coeff;
	PR_VEL(:,5) = PR_VEL(:,5)*Coeff;
	PR_VEL(:,6) = PR_VEL(:,6)*Coeff;
	PR_VEL(:,7) = PR_VEL(:,7)*Coeff;
	PR_VEL(:,8) = PR_VEL(:,8)*Coeff;
	PR_VEL(:,9) = PR_VEL(:,9)*Coeff;
	PR_VEL(:,10) = PR_VEL(:,10)*Coeff;
	PR_VEL(:,11) = PR_VEL(:,11)*Coeff;
	

	
	for m=5:11
        PR_VEL(:,m) = PR_VEL(:,m);
	end
        
	for m = 1:6
        subplot(2,3,m); p = plot(PR_VEL(:,1),PR_VEL(:,m+1),'.'); suc = PRespTitle(m,p,2);xlabel('time (s)');ylabel('uM m-2 s-1');
	end
	pause;
	
	for m = 5:10
        subplot(2,3,m-4); p = plot(PR_VEL(:,1),PR_VEL(:,m+1),'.'); suc = PRespTitle(m,p,2);xlabel('time (s)');ylabel('uM m-2 s-1');
	end
	pause;
	
	% Get the fluxes for reactions in photosynthesis pathway
	
	PS_VEL = PS_VEL*Coeff;
	PS_VEL(1,:)=PS_VEL(1,:)/(Coeff);
	
	for m = 1:6
	subplot(2,3,m); p = plot(PS_VEL(1,:), PS_VEL(1+m,:),'.'); suc = PSTitle(m,p,2);
	end
	pause;
	
	for m = 7:12
	subplot(2,3,m-6); p = plot(PS_VEL(1,:), PS_VEL(1+m,:),'.'); suc = PSTitle(m,p,2);
	end
	pause;
	
	for m = 12:17
	subplot(2,3,m-11); p = plot(PS_VEL(1,:), PS_VEL(1+m,:),'.'); suc = PSTitle(m,p,2);
	end
	pause;
	
	% The following plot the net carbon assimilation rate based on PS_VEL and PR_VEL.
	
	CarbonRate = PS_VEL(2,:);
	CO2Release = PR_VEL(:,9);
	NetCarbon = CarbonRate' - CO2Release;
	subplot(1,1,1);
	plot(PS_VEL(1,:),NetCarbon,'.');
	
	title('The net A from PS_PROUT (NOT USED IN DynaPS)');
	xlabel('time(second)');
	ylabel('umol m-2 s-1');
	pause; 
	
	PSVEL = PS_VEL';
	PRVEL = PR_VEL';
	
	global CO2A; 
	
	% CO2A(:,4) = NetCarbon;
	
	PS_VEL(2:18,:) = PS_VEL(2:18,:)/Coeff; 
	PR_VEL(:,2:11) = PR_VEL(:,2:11)/Coeff; 


Success = 1; 

end