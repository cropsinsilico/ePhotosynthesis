function out = FIBF_Dict(cond, start)

global FIBF_Pool;
global FIBF_AUX;

pool_values = [
FIBF_Pool(1), ...
];

out = dictionary();

out_MOD = import_pt_dict("FIBF", "MOD");
out_COND = import_pt_dict("FIBF", "COND");
out_POOL = dictionary(FIBF_NAMES("POOL"), pool_values);

out_MOD('FIBF2FI_PQa') = FIBF_AUX(2);

out('COND') = out_COND;
out('POOL') = out_POOL;
out('MOD') = out_MOD;

end
