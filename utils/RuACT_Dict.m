function out = RuACT_Dict(cond, start)

out = dictionary();

global RuACT_Pool;
global RuACT_RC;

out_MOD = import_pt_dict("RuACT", "MOD");
out_COND = import_pt_dict("RuACT", "COND");
out_POOL = dictionary(RuACT_NAMES("POOL"), RuACT_Pool');
out_RC = dictionary(RuACT_NAMES("RC"), RuACT_RC');

out('COND') = out_COND;
out('MOD') = out_MOD;
out('POOL') = out_POOL;
out('RC') = out_RC;

out_VEL = import_pt_dict("RuACT", "VEL");
if (numEntries(out_VEL) > 0)
    out('VEL') = out_VEL;
end

end
