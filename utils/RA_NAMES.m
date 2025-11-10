function out = RA_NAMES(pt)
global RA_COND_START;
RA_COND_START = 92;
if (pt == "CHILDREN")
out = ["EPS", "RuACT"];
elseif (pt == "PARAM_TYPES")
out = [];
elseif (pt == "MOD")
out = [
];
elseif (pt == "COND")
out = [
];
end
end
