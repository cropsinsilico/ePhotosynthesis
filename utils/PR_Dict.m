function out = PR_Dict(cond, start)

out = dictionary();

out_COND = import_pt_dict("PR", "COND");
PS_COND = import_pt_dict("PS", "COND");
out_COND('RuBP') = PS_COND('RuBP');
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
