function out = export_cond_dict(mod, cond)
global export_mod_enabled;
if ~export_mod_enabled
    out = dictionary();
    return;
end
names = eval(mod + '_NAMES(''COND'');');
eval("global " + mod + "_COND_START;");
eval("start = " + mod + "_COND_START;");
values = reshape(cond(start:(start + length(names) - 1)), 1, []);
out = dictionary(names, values);
export_pt_dict(mod, "COND", out);
end
