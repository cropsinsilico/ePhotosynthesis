function out = PS_Dict(cond, start)

out = dictionary();

out_COND = import_pt_dict("PS", "COND");
out_COND('ADPG') = 0;
out_COND('_Pi') = 0;
out_COND('_ADP') = 0;
out_COND('_v1') = 0;

out_MOD = import_pt_dict("PS", "MOD");
out_MOD('_NADPH') = out_COND('NADPH');
out_MOD('RegFactor') = 1;
out_MOD('SC') = 1;
out_MOD('SC1') = 1;
out_MOD('STOM1') = 1;
out_MOD('STOM2') = 1;
out_MOD = remove(out_MOD, 'KE1Ratio');
out_MOD = remove(out_MOD, 'KE2Ratio');

out_COND = remove(out_COND, 'ADPG');
out_COND = remove(out_COND, 'NADPH');
out_COND = remove(out_COND, 'CO2');
out_COND = remove(out_COND, 'O2');
out_MOD = remove(out_MOD, 'KE57'); % Only used for C3 in C++

out('COND') = out_COND;
out('MOD') = out_MOD;

out_VEL = import_pt_dict("PS", "VEL");
if (numEntries(out_VEL) > 0)
    out('VEL') = out_VEL;
end

end
