function out = FIBF_NAMES(pt)
global FIBF_COND_START;
FIBF_COND_START = 52;
if (pt == "CHILDREN")
out = ["BF", "FI"];
elseif (pt == "PARAM_TYPES")
out = ["pool"];
elseif (pt == "MOD")
out = [
"ChlPSI", ...
"ChlT", ...
"ChlT2", ...
"FIBF2FI_PQ", ...
"FIBF2FI_PQa", ...
];
elseif (pt == "COND")
out = [
"kd", ...
];
elseif (pt == "POOL")
out = [
"PQT", ...
];
end
end