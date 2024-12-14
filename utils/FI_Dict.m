function out = FI_Dict(cond, start)

out = dictionary();

global FI_Pool;
global FI_RC;

out_MOD = import_pt_dict("FI", "MOD");
out_COND = import_pt_dict("FI", "COND");
out_POOL = dictionary(FI_NAMES("POOL"), FI_Pool);
out_RC = dictionary(FI_NAMES("RC"), FI_RC');

out('COND') = out_COND;
out('POOL') = out_POOL;
out('MOD') = out_MOD;
out('RC') = out_RC;

out_VEL = import_pt_dict("FI", "VEL");
if (numEntries(out_VEL) > 0)
    out('VEL') = out_VEL;
end

end
