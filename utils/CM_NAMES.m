function out = CM_NAMES(pt)
global CM_COND_START;
CM_COND_START = 36;
if (pt == "CHILDREN")
out = ["PS_PR", "SUCS"];
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
