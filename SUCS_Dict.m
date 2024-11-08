function out = SUCS_Dict(CM_Con)

out = dictionary();

global SUCS_Pool;

out_COND = import_pt_dict("SUCS", "COND");
out_POOL = dictionary(SUCS_NAMES("POOL"), SUCS_Pool);

out_MOD = import_pt_dict("SUCS", "MOD");
out_MOD = remove(out_MOD, 'KI583');
out_MOD('SUCS2PS_Pic') = out_MOD('PsPEXT');
out_MOD = remove(out_MOD, 'PsPEXT');
out_MOD = remove(out_MOD, 'SC');
out_MOD = remove(out_MOD, 'SC1');
out_MOD = remove(out_MOD, 'KE5Ratio');

out_COND = remove(out_COND, 'UTPc');
out_COND = remove(out_COND, 'ATPc');
out_COND = remove(out_COND, 'ADPc');
out_COND = remove(out_COND, 'OPOPc');

out('COND') = out_COND;
out('POOL') = out_POOL;
out('MOD') = out_MOD;

out_VEL = import_pt_dict("SUCS", "VEL");
if (numEntries(out_VEL) > 0)
    out('VEL') = out_VEL;
end

end
