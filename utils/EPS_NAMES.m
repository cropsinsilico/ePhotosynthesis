function out = EPS_NAMES(pt)
global EPS_COND_START;
EPS_COND_START = 87;
if (pt == "CHILDREN")
out = ["FIBF", "PS_PR"];
elseif (pt == "PARAM_TYPES")
out = [];
elseif (pt == "MOD")
out = [];
elseif (pt == "COND")
out = [];
elseif (pt == "VEL")
out = [];
end
end
