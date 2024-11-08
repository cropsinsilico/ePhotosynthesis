function out = export_mod_dict(mod, cond, prefix)
global export_mod_enabled;
if ~export_mod_enabled
    out = dictionary();
    return;
end
% TODO: VEL or just on import side?
export_cond_dict(mod, cond);
names = eval(mod + '_NAMES(''MOD'');');
values = zeros(size(names));
for i = 1:length(names)
    iname = names{i};
    ival = get_global_var(iname);
    if isempty(ival)
        try
            ival = evalin('caller', iname);
        catch ME
            if (~(strcmp(ME.identifier,'MATLAB:UndefinedFunction') | ...
                    strcmp(ME.identifier, 'MATLAB:m_illegal_character')))
                rethrow(ME);
            end
        end
    end
    if isempty(ival) && strlength(prefix) > 0
        ival = get_global_var(prefix + iname);
        if isempty(ival)
            try
                ival = evalin('caller', prefix + iname);
            catch ME
                if (~(strcmp(ME.identifier,'MATLAB:UndefinedFunction') | ...
                        strcmp(ME.identifier, 'MATLAB:m_illegal_character')))
                    rethrow(ME);
                end
            end
        end
    end
    if ~isempty(ival)
        values(i) = ival;
    end
end
out = dictionary(names, values);
eval("global " + mod + "_MOD_DICT;");
eval(mod + "_MOD_DICT = out;");
end
