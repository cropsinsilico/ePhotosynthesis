function out = trDynaPS_NAMES(pt)
global trDynaPS_COND_START;
trDynaPS_COND_START = 120;
if (pt == "CHILDREN")
out = ["DynaPS", "RROEA"];
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
