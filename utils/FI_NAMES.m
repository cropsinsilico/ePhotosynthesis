function out = FI_NAMES(pt)
global FI_COND_START;
FI_COND_START = 1;
if (pt == "CHILDREN")
out = [];
elseif (pt == "PARAM_TYPES")
out = ["POOL", "RC"];
elseif (pt == "MOD")
out = [
"cpsii", ...
];
elseif (pt == "COND")
out = [
"A", ... % 	FI_Con	(	1	)			;	% 	The concentration of excitons in the peripheral antenna
"U", ... % 	FI_Con	(	2	)			;	%	The concentration fo excitons in the core antenna
"P680ePheo", ... %    FI_Con  (   3   );          %QF add
"P680pPheon", ... % 	FI_Con	(	4	)			;	%	The concentration for the P680+ Pheo-
"P680pPheo", ... % 	FI_Con	(	5	)			;	% 	The concentration of P680+ Pheo
"P680Pheon", ... % 	FI_Con	(	6	)			;	%	The concentration of P680Pheo-
"Yz", ... % 	FI_Con	(	7	)			;	%	The concentration of reduced tyrosine
"S1T", ... % 	FI_Con	(	8	)			;	% 	The concentration of S1 in combination with reduced tyrosine 
"S2T", ... % 	FI_Con	(	9	)			;	%	The concentration of S2 in combination with reduced tyrosine
"S3T", ... % 	FI_Con	(	10	)			;	% 	The concentration of S3 in combination with reduced tyrosine
"S0T", ... % 	FI_Con	(	11	)			;	%	The concentration of S0 in combination with reduced tyrosine 
"S1Tp", ... % 	FI_Con	(	12	)		;	% 	The concentration of S1 in combination with oxidized tyrosine
"S2Tp", ... % 	FI_Con	(	13	)		;	% 	The concentration of S2 in combination with oxidized tyrosine
"S3Tp", ... % 	FI_Con	(	14	)		;	%	The concentration of S3 in combination with oxidized tyrosine
"S0Tp", ... % 	FI_Con	(	15	)		;	% 	The concentration of S0 in combination with oxidized tyrosine
"QAQB", ... % 	FI_Con	(	16	)		;	% 	The concentration of [QAQB]
"QAnQB", ... % 	FI_Con	(	17	)		;	% 	The concentration of [QA-QB];
"QAQBn", ... % 	FI_Con	(	18	)		;	%	The concentration of [QAQB-]
"QAnQBn", ... % 	FI_Con	(	19	)		;	% 	The concentration of [QA-QB-];
"QAQB2n", ... % 	FI_Con	(	20	)		;	%	The concentration of [QAQB2-]
"QAnQB2n", ... % 	FI_Con	(	21	)		;	% 	The concentration of [QA-QB2-];
"PQn", ... % 	FI_Con	(	22	)			;	%	The concentration of reduced PQ, i.e. PQH2;				
];
elseif (pt == "POOL")
out = [
"QBt", ...
"PQT", ...
];
elseif (pt == "RC")
out = [
"kA_d", ... % 	FI_RC	(	1	)				;	%	The rate constant of heat dissipation from peripheral antenna	Lazar (1999), 0.25~1 *10^(9)			
"kA_f", ... % 	FI_RC	(	2	)				;	%	The rate constant of fluorescence emission from peripheral antenna	Lazar 1999, with a lifetime of 5 ns at closed reaction center			
"kA_U", ... % 	FI_RC	(	3	)				;	%	The rate constant of exciton transfer from periphral antenna to core antenna	Reference needed, a guess			
"kU_A", ... % 	FI_RC	(	4	)				;	%	The rate constant of exciton transfer from core antenna to peripheral antenna	Reference needed, a guess			
"kU_d", ... % 	FI_RC	(	5	)				;	%	The rate constant of  heat emission from core antenna				
"kU_f", ... % 	FI_RC	(	6	)				;	%	The rate constant of fluorescence emission from core antenna				
"k1", ... % 	FI_RC	(	7	)		        ;	%	The rate constant of primary charge separation for open reaction center				
"k_r1", ... % 	FI_RC	(	8	)				;	%	The rate constant of charge recombination for open reactoin center				
"kz", ... % 	FI_RC	(	9	)				    ;	%	The rate constant of the Tyrosine oxidation	Lazar (1999); 3.8~50 * 10^6 			
"k12", ... % 	FI_RC	(	10	)				    ;	%	The rate constant of the S1 to S2 transition	Lazar (1999); 0.667~33.3 * 10^3			
"k23", ... % 	FI_RC	(	11	)				    ;	%	The rate constant of the S2 to S3 transition	Lazar (1999); 0.667~33.3 * 10^3			
"k30", ... % 	FI_RC	(	12	)				    ;	%	The rate constant of the S3 to S0 transition	Lazar (1999); 0.667~33.3 * 10^3			
"k01", ... % 	FI_RC	(	13	)				    ;	%	The rate constant of the S0 to S1 transition	Lazar (1999); 0.667~33.3 * 10^3			
"k2", ... % 	FI_RC	(	14	)				    ;	%	The rate constant of the QA reduction by Pheo-	Lazar (1999); 2~2.3 * 10^9			
"kAB1", ... % 	FI_RC	(	15	)				;	%	The rate constant of QAQB-->QAQB-	Lazar (1999); 2.5~5 * 10^3			
"kBA1", ... % 	FI_RC	(	16	)				;	%	The rate constant of the QAQB- -->QAQB	Lazar (1999)			
"kAB2", ... % 	FI_RC	(	17	)				;	%	The rate constant of the QAQB- --> QAQB2-	Lazar (1999); 1.25~3.33 * 10^3			
"kBA2", ... % 	FI_RC	(	18	)				;	%	The rate constant of the QAQB2- --> QAQB- 	Lazar (1999), or same as kAB2 depend on the equilibium constant			
"k3", ... % 	FI_RC	(	19	)				    ;	%	The rate constant of the exchange of PQ and QBH2	Lazar (1999),0.12~1 for the fast PQ pool,  or 3~8 for the slow recycling PQ pool			
"k_r3", ... % 	FI_RC	(	20	)				;	%	The rate constant of the exchange of QB and PQH2	Lazar (1999), since the equilibrium constant is 1 (205 in Lazar, 1999) 			
"k_pq_oxy", ... % 	FI_RC	(	21	)				;	%	The rate constant of the PQH2 oxidation	Lazar (1999),50~500			
];
elseif (pt == "VEL")
out = [
"vA_d", ... %	vA_d	The rate of heat dissipation from peripheral antenna
"vA_f", ... %	vA_f	The rate of fluorescence emission from peripheral antenna
"vA_U", ... %	vA_U	The rate of exciton transfer from peripheral antenna to core antenna in open reaction center
"vU_A", ... %	vU_A	The rate of exciton transfer from core antenna to perpheral antenna in open center
"vU_f", ... %	vU_f	The rate of fluorescence emission from core antenna
"vU_d", ... %	vU_d	The rate of heat dissipation from core antenna
"v1", ... %	v1	The rate of primary charge separation
"v_r1", ... %	v_r1	The rate of charge recombination
"vS1_S2", ... %	vS1_S2	The rate of transition from S1 to S2 
"vS2_S3", ... %	vS2_S3	The rate of transition from S2 to S3
"vS3_S0", ... %	vS3_S0	The rate of transition from S3 to S0
"vS0_S1", ... %	vS0_S1	The rate of transition from S0 to S1
"vz_1", ... %	vz_1	The rate of P680p reduction
"v1z_1", ... %	v1z_1	The rate of oxidation of S1T by P680pPheon
"v2z_1", ... %	v2z_1	The rate of oxidation of S2T  by P680pPheon
"v3z_1", ... %	v3z_1	The rate of oxidation of S3T  by P680pPheon
"v0z_1", ... %	v0z_1	The rate of oxidation of S0T  by P680pPheon
"vz_2", ... %	vz_2	The rate of P680pPheon reduction
"v1z_2", ... %	v1z_2	The rate of oxidation of S1T by P680pPheo
"v2z_2", ... %	v2z_2	The rate of oxidation of S2T  by P680pPheo
"v3z_2", ... %	v3z_2	The rate of oxidation of S3T  by P680pPheo
"v0z_2", ... %	v0z_2	The rate of oxidation of S0T  by P680pPheo
"v1z", ... %	v1z	
"v2z", ... %	v2z	
"v3z", ... %	v3z	
"v0z", ... %	v0z	
"vAB1", ... %	vAB1	The rate of electron transfer from QA- to QB
"vBA1", ... %	vBA1	The rate of electron transfer from QB- to QA
"vAB2", ... %	vAB2	The rate of electron transfer from QA- to QB-
"vBA2", ... %	vBA2	The rate of electron transfer from QB2- TO QA
"v3", ... %	v3	The rate of exchange of QAQBH2 with PQ
"v_r3", ... %	v_r3	The rate of exchange of QAQB with PQH2
"v3_n", ... %	v3_n	The rate of exchange of QAnQBH2 with PQ
"v_r3_n", ... %	v_r3_n	The rate of exchange of QAnQB with PQH2
"v_pq_ox", ... %	v_pq_ox	The rate of PQH2 oxidation
"Ic", ... %	Ic	The incident light on the core antenna
"Ia", ... %	Ia	The incident light on the peripheral antenna
"v2_1", ... %	v2_1	The rate of P680pPheon oxidation
"v2_2", ... %	v2_1	The rate of P680pPheon oxidation
"v2_00_1", ... %	v2_00_1	The rate of reduction of QAQB by P680pPheon
"v2_01_1", ... %	v2_01_1	The rate of reduction of QAQBn by P680pPheon
"v2_02_1", ... %	v2_02_1	The rate of reduction of QAQB2n by P680pPheon
"v2_00_2", ... %	v2_00_2	The rate of reduction of QAQB by P680Pheon
"v2_01_2", ... %	v2_01_2	The rate of reduction of QAQBn by P680Pheon
"v2_02_2", ... %	v2_02_2	The rate of reduction of QAQB2n by P680Pheon
"vr2_00_1", ... %	vr2_00_1	The reverse reaction of The rate of reduction of QAQB by P680pPheon
"vr2_01_1", ... %	vr2_01_1	The reverse reaction of The rate of reduction of QAQBn by P680pPheon
"vr2_02_1", ... %	vr2_02_1	The reverse reaction of The rate of reduction of QAQB2n by P680pPheon
"vr2_1", ... %	vr2_1	
"vr2_00_2", ... %	vr2_00_2	The reverse reaction of The rate of reduction of QAQB by P680Pheon
"vr2_01_2", ... %	vr2_01_2	The reverse reaction of The rate of reduction of QAQBn by P680Pheon
"vr2_02_2", ... %	vr2_02_2	The reverse reaction of The rate of reduction of QAQB2n by P680Pheon
"vr2_2", ... %	vr2_2	
"vP680qU", ... %	vr2_2	
"vP680qA", ... %	vr2_2	
"vU_P680", ... 
"vP680_d", ... 
"vP680_f", ... 
];
end
end
