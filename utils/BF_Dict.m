function out = BF_Dict(cond, start)

out = dictionary();

global BF_Pool;
global BF_RC;

out_MOD = import_pt_dict("BF", "MOD");
out_COND = import_pt_dict("BF", "COND");
out_POOL = dictionary(BF_NAMES("POOL"), BF_Pool');
out_RC = dictionary(BF_NAMES("RC"), BF_RC');

out_MOD('_Pi') = out_COND('Pi');
out_MOD('PMODTEM') = 1.0;
out_MOD('RT') = 298 * 8.314;
out_MOD('DeltaGo') = 7.3 * 4184;
out_COND = remove(out_COND, "Pi");

out('COND') = out_COND;
out('POOL') = out_POOL;
out('MOD') = out_MOD;
out('RC') = out_RC;

out_VEL = import_pt_dict("BF", "VEL");
if (numEntries(out_VEL) > 0)
    out('VEL') = out_VEL;
end

end
