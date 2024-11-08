function out = export_pt_dict(mod, pt, x)
global export_mod_enabled;
if ~export_mod_enabled
    out = dictionary();
    return;
end
out = x;
eval("global " + mod + "_" + pt + "_DICT;");
eval(mod + "_" + pt + "_DICT = out;");
end
