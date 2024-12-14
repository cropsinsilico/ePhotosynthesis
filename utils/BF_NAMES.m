function out = BF_NAMES(pt)
global BF_COND_START;
BF_COND_START = 1;
if (pt == "CHILDREN")
out = [];
elseif (pt == "PARAM_TYPES")
out = ["POOL", "RC"];
elseif (pt == "MOD")
out = [
"cATPsyn", ...
"CPSi", ...
"cNADPHsyn", ...
"EPS_ATP_Rate", ...
"_Pi", ...
"PMODTEM", ...
"RT", ...
"DeltaGo", ...
];
elseif (pt == "COND")
out = [
"ISPHr", ... %	BF_con	(	1	)			;	%	The reduced ion sulfer protein (ISPH)
"cytc1", ... %	BF_con	(	2	)			;	%	The oxidized state of cytc1
"ISPo", ... %	BF_con	(	3	)			;	%	The oxidized ion sulfer protein (ISP)
"ISPoQH2", ... %	BF_con	(	4	)			;	%	The complex of oxidized ion sulfer protein and reduced quinone
"QHsemi", ... %	BF_con	(	5	)			;	%	Semiquinone
"cytbL", ... %	BF_con	(	6	)		;	%	The oxidized cytbL
"Qi", ... %	BF_con	(	7	)			;	%	The binding Quinone
"Q", ... %	BF_con	(	8	)			;	%	Quinone
"cytbH", ... %	BF_con	(	9	)		;	%	The oxidized form of cytbH
"Qn", ... %	BF_con	(	10	)			;	%	Q-
"Qr", ... %	BF_con	(	11	)			;	%	Q2-
"QH2", ... %	BF_con	(	12	)			;	%	QH2
"cytc2", ... %	BF_con	(	13	)		;	%	oxidized cytc2
"P700", ... %	BF_con	(	14	)		;	%	The reduced state of P700, including both P700 and excited P700
"ADP", ... %	BF_con	(	15	)			;	%	ADP in stroma
"Pi", ... %	BF_con	(	16	)			;	%	Phosphate in stroma
"ATP", ... %	BF_con	(	17	)			;	%	ATP in stroma
"Ks", ... %	BF_con	(	18	)			;	%	K ions in stroma
"Mgs", ... %	BF_con	(	19	)			;	%	Mg ions in stroma
"Cls", ... %	BF_con	(	20	)			;	%	Cl ions in stroma
"Aip", ... %	BF_con	(	21	)			;	%	The number of photons in peripheral antenna
"U", ... %	BF_con	(	22	)			;	%	The number of photons in core antenna
"An", ... %	BF_con	(	23	)			;	%	The reduced electron acceptor in PSI
"Fdn", ... %	BF_con	(	24	)			;	%	The reduced ferrodoxin
"BFHs", ... %	BF_con	(	25  )		    ;	%	The proton concentration in stroma, put in unit: mol l-1
"BFHl", ... %	BF_con	(	26	)		    ;	%	The proton concentration in lumen; put in unit: mol l-1
"PHs", ... %	BF_con	(	27	)			;	%	The PH value of stroma
"PHl", ... %	BF_con	(	28	)			;	%	The PH value of lumen
"NADPH", ... %	BF_con	(	29	)		    ;	%	The NADPH concentratin in stroma, Unit: mmol l-1
];
elseif (pt == "VEL")
out = [
"Vbf1", ... % BF_Vel	(	1	)	=	Vbf1	;
"Vbf2", ... % BF_Vel	(	2	)	=	Vbf2	;
"Vbf3", ... % BF_Vel	(	3	)	=	Vbf3	;
"Vbf4", ... % BF_Vel	(	4	)	=	Vbf4	;
"Vbf5", ... % BF_Vel	(	5	)	=	Vbf5	;
"Vbf6", ... % BF_Vel	(	6	)	=	Vbf6	;
"Vbf7", ... % BF_Vel	(	7	)	=	Vbf7	;
"Vbf8", ... % BF_Vel	(	8	)	=	Vbf8	;
"Vbf9", ... % BF_Vel	(	9	)	=	Vbf9	;
"Vbf10", ... % BF_Vel	(	10	)	=	Vbf10	;
"Vbf11", ... % BF_Vel	(	11	)	=	Vbf11	;
"Vqi", ... % BF_Vel	(	12	)	=	Vqi	;
"Vipc", ... % BF_Vel	(	13	)	=	Vipc	;
"Vicp", ... % BF_Vel	(	14	)	=	Vicp	;
"Vinc", ... % BF_Vel	(	15	)	=	Vinc	;
"Vinp", ... % BF_Vel	(	16	)	=	Vinp;
"Vdp", ... % BF_Vel	(	17	)	=	Vdp	;
"Vdc", ... % BF_Vel	(	18	)	=	Vdc	;
"Vfp", ... % BF_Vel	(	19	)	=	Vfp	;
"Vfc", ... % BF_Vel	(	20	)	=	Vfc	;
"Vsfd", ... % BF_Vel	(	21	)	=	Vsfd	;
"VsATP", ... % BF_Vel	(	22	)	=	VsATP	;
"VgPQH2", ... % BF_Vel	(	23	)	=	VgPQH2	;
"Vbf12", ... % BF_Vel	(	24	)	=	Vbf12	;
"Vbf13", ... % BF_Vel	(	25	)	=	Vbf13	;
"Vbf14", ... % BF_Vel	(	26	)	=	Vbf14	;
"Vbf15", ... % BF_Vel	(	27	)	=	Vbf15	;
"Vbf16", ... % BF_Vel	(	28	)	=	Vbf16	;
"vbfn2", ... % BF_Vel	(	29	)	=	vbfn2	;
"VsNADPH", ... % BF_Vel	(	30	)	=	VsNADPH	;
"vcet", ... % BF_Vel	(	31	)	=	vcet ;
];
elseif (pt == "POOL")
out = [
"Tcyt", ... %	BF_Pool	(	1	)			;	%	The total amount of cytbH or cytbL; Unit: micromole m-2 leaf area				
"Tcytc2", ... %	BF_Pool	(	2	)			;	%	The total amount of cytc; Unit: micromole m-2 leaf area				
"TK", ... %	BF_Pool	(	3	)			;	%	The total concentration of K in both stroma and lumen. Unit: mmol l-1. In this model, it was assumed that the total concentration of K, and Mg and Cl as well, is constant. 				
"TMg", ... %	BF_Pool	(	4	)			;	%	The total concentration of Mg in both stroma and lumen. Unit: mmol l-1. In this model, it was assumed that the total concentration of Mg, and K and Cl as well, is constant. 				
"TCl", ... %	BF_Pool	(	5	)			;	%	The total concentration of Cl in both stroma and lumen. Unit: mmol l-1. In this model, it was assumed that the total concentration of Cl in both stroma and lumen is constant. 				
"TFd", ... %	BF_Pool	(	6	)			;	%	The total concentration of Ferrodoxin				
"TA", ... %	BF_Pool	(	7	)			;	%	The total concentration of the primary electron acceptor of PSI; Unit: micromole m-2 leaf area				
"TQ", ... %	BF_Pool	(	8	)			;	%	The total concentration of plastoquinone in thylakoid membrane. ; Unit: micromole m-2 leaf area				
"BFTs", ... %	BF_Pool	(	9	)			;	%	The total concentration of buffer in stroma; unit: mmol per liter				
"BFTl", ... %	BF_Pool	(	10	)			;	%	The total concentration of buffer in lumen; unit: mmol per liter	
"P700T", ... %	BF_Pool	(	11	)			;	%	The total concentration of P700; unit: micromole m-2 leaf area				
"NADPHT", ... %   BF_Pool (12)                ;   %   The total concentration of NADPH, Unit: mmol l-1.
];
elseif (pt == "RC")
out = [
"K1", ... %	BF_RC	(	1	)			;	%	The rate constant for formation of ISP.QH2 complex; unit:  per second 
"K2", ... %	BF_RC	(	2	)			;	%	The rate constant for ISP.QH2-->QH(semi) + ISPH(red) ; unit:  per second	
"K3", ... %	BF_RC	(	3	)			;	%	The rate constant for QH. + cytbL --> Q + cytbL- + H+	Unit: s-1
"K4", ... %	BF_RC	(	4	)			;	%	The rate constant for cytbL- + cytbH --> cytbL + cytbH-	Unit: s-1
"K5", ... %	BF_RC	(	5	)			;	%	The rate constant for CytbH- + Q --> cytbH + Q-	Unit: s-1
"K6", ... %	BF_RC	(	6	)			;	%	The rate constant  for CytbH- + Q- --> cytbH + Q2-	Unit: s-1
"K7", ... %	BF_RC	(	7	)			;	%	The rate constant for Q binding to Qi site; which assumed half time as 200 us, following Croft's website	Unit: s-1
"K8", ... %	BF_RC	(	8	)			;	%	The rate constant for ISPH + CytC1 --> ISPH(ox) + CytC1+	Unit: s-1
"K9", ... %	BF_RC	(	9	)			;	%	The rate constant for the electron transport from cytc1 to cytc2	Unit: s-1
"K10", ... %	BF_RC	(	10	)			;	%	The rate constant for the electron transport from cytc2 to P700	Unit: s-1
"Vmax11", ... % BF_RC	(	11	)			;	%	The maximum rate of ATP synthesis	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1
"Kqi", ... %	BF_RC	(	12	)			;	%	The rate constant for uptake of two protons from the stroma to Q2-	s-1
"PK", ... %	BF_RC	(	13	)			;	%	The permeability constant for K	Unit: cm s-1
"PMg", ... %	BF_RC	(	14	)			;	%	The permeability constant for Mg	Unit: cm s-1
"PCl", ... %	BF_RC	(	15	)			;	%	The permeability constant for Cl	Unit: cm s-1
"Kau", ... %	BF_RC	(	16	)			;	%	The rate constant for exciton transfer from perpheral antenna to core antenna, see FI	Unit: s-1
"Kua", ... %	BF_RC	(	17	)			;	%	The rate constant for exciton transfer from core antenna to peripheral antenna, SEE FI	Unit: s-1
"Kf", ... %	BF_RC	(	18	)			;	%	The rate constant for fluorescence emission, see the note in FI	Unit: s-1
"Kd", ... %	BF_RC	(	19	)			;	%	The rate constant for heat dissipation; see the note for FI	Unit: s-1
"KE8", ... %	BF_RC	(	20	)			;	%	ISPHr + cytc1 --> ISPHox + cytc1-	Unit: s-1
"KE9", ... %	BF_RC	(	21	)			;	%	cytc1- + cytc2 --> cytc1 + cytc2-	Unit: s-1
"K15", ... %	BF_RC	(	22	)			;	%	The rate constant for primary charge separation in PSI	Unit: s-1
"K16", ... %	BF_RC	(	23	)			;	%	The rate constant for electron tranfer from electron acceptor of PSI to Fd	Unit: s-1
"MemCap", ... %	BF_RC(	24	)			;	%	The membrane capacity	
"RVA", ... %	BF_RC	(	25	)			;	%	The ratio of lumen volume to thylakoid membrane area	
"KBs", ... %	BF_RC	(	26	)			;	%	The buffer equilibrium constant in stroma	
"KBl", ... %	BF_RC	(	27	)			;	%	The buffer equilibrium constant in lumen	
"KM1ATP", ... %	BF_RC	(	28	)			;	%	The michaelis menton constant for ATP for ATP synthesis
"KM1ADP", ... %	BF_RC	(	29	)			;	%	The michaelis menton constant for ATP for ADP synthesis
"KM1PI", ... %	BF_RC	(	30	)			;	%	The michaelis menton constant for ATP for PI synthesis
"KM2NADP", ... %	BF_RC	(	31	)			;	%	The michaelis menten constant for NADP	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1
"KM2NADPH", ... %   BF_RC(	32	)			;	%	The michaelis menten constant for NADPH	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1
"V2M", ... %	BF_RC	(	33	)			;	%	The maximum rate of NADPH formation	Unit: mmol l-1 s-1; The unit for the reactions occurrs in stroma is mmol l-1 s-1
"KE2", ... %	BF_RC	(	34	)			;	%	Equilibrium constatn	
];
end
end
