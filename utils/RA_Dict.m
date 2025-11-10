function out = RA_Dict(cond, start)

out = dictionary();

out_MOD = import_pt_dict("RA", "MOD");
out_COND = import_pt_dict("RA", "COND");

out('COND') = out_COND;
out('MOD') = out_MOD;

end
