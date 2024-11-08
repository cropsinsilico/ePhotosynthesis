function out = PR_Dict(CM_Con)

out = dictionary();

out_COND = import_pt_dict("PR", "COND");
out_COND('RuBP') = CM_Con(1);
out_COND('_v131') = 0;

out_MOD = import_pt_dict("PR", "MOD");
out_MOD('CE') = 1;

out('COND') = out_COND;
out('MOD') = out_MOD;

out_VEL = import_pt_dict("PR", "VEL");
if (numEntries(out_VEL) > 0)
    out('VEL') = out_VEL;
end

end
