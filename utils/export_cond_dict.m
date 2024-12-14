function out = export_cond_dict(mod, cond, start)
arguments
    mod (1,1) string
    cond (:,:) double
    start (1,1) double = 1
end
global export_mod_enabled;
if ~export_mod_enabled
    out = dictionary();
    return;
end
names = eval(mod + '_NAMES(''COND'');');
values = reshape(cond(start:(start + length(names) - 1)), 1, []);
if isempty(names) && isempty(values)
    out = dictionary();
else
    out = dictionary(names, values);
end
export_pt_dict(mod, "COND", out);
end
