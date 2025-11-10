function out = DynaPS_NAMES(pt)
global DynaPS_COND_START;
DynaPS_COND_START = 96;
if (pt == "CHILDREN")
out = ["RA", "XanCycle"];
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
