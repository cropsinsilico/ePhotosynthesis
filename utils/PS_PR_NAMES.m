function out = PS_PR_NAMES(pt)
global PS_PR_COND_START;
PS_PR_COND_START = 1;
if (pt == "CHILDREN")
out = ["PS", "PR"];
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
