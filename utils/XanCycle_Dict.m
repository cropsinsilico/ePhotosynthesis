function out = XanCycle_Dict(cond, start)

global XanCycle_kva;
global XanCycle_kaz;
global XanCycle_kza;
global XanCycle_kav;

out = dictionary();

out_MOD = import_pt_dict("XanCycle", "MOD");
out_COND = import_pt_dict("XanCycle", "COND");
out_MOD('Vx_') = 160;
out_MOD('Ax_') = 10;
out_MOD('Zx_') = 5;
out_MOD('ABA_') = 1;

out_MOD('kva') = XanCycle_kva;
out_MOD('kaz') = XanCycle_kaz;
out_MOD('kza') = XanCycle_kza;
out_MOD('kav') = XanCycle_kav;

out('COND') = out_COND;
out('MOD') = out_MOD;

out_VEL = import_pt_dict("XanCycle", "VEL");
if (numEntries(out_VEL) > 0)
    out('VEL') = out_VEL;
end

end
