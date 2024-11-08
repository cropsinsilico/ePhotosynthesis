function out = CM_Dict(CM_Con)
out = dictionary();
out('PR') = PR_Dict(CM_Con);
out('PS') = PS_Dict(CM_Con);
out('SUCS') = SUCS_Dict(CM_Con);
out('ALL') = VARS_Dict();
end
