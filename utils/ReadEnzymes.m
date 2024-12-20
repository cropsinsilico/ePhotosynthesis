function out = ReadEnzymes(name)
order = [
"V1", ...
"V2", ...
"V3", ...
"V5", ...
"V6", ...
"V7", ...
"V8", ...
"V9", ...
"V10", ...
"V13", ...
"V23", ...
"V16", ...
"V112", ...
"V113", ...
"V121", ...
"V122", ...
"V123", ...
"V124", ...
"V131", ...
"V51", ...
"V52", ...
"V55", ...
"V56", ...
"V57", ...
"V58", ...
"V59", ...
"Jmax", ...
"K1", ...
"K2", ...
"K3", ...
"K4", ...
"K5", ...
"K6", ...
"K7", ...
"K8", ...
"K9", ...
"K10", ...
"Vmax11", ...
"Kau", ...
"Kua", ...
"Kf", ...
"Kd ", ...
"K15", ...
"K16", ...
"V2M ", ...
"kA_d", ...
"kA_f", ...
"kA_U", ...
"kU_A", ...
"kU_d", ...
"kU_f", ...
"k1", ...
"k_r1", ...
"kz", ...
"k12", ...
"k23", ...
"k30", ...
"k01", ...
"k2", ...
"kAB1", ...
"kBA1", ...
"kAB2", ...
"kBA2", ...
"k3", ...
"k_r3", ...
"k_pq_oxy", ...
];
InputF = importdata(name);
keys = InputF.textdata(:, 1);
out = ReorderData(order, keys, InputF.data(:, 1), 0.0);
end
