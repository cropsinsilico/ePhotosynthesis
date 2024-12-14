function out = export_mod_dict(t, mod, cond, prefix, parent, start)
arguments
    t (1,1) double
    mod (1,1) string
    cond (:,:) double
    prefix (1,1) string = ""
    parent (1,1) string = ""
    start (1,1) double = 1
end
global export_mod;
global export_mod_out;
global export_mod_enabled;
if isempty(export_mod_enabled) || ~export_mod_enabled
    out = dictionary();
    return;
end
% TODO: VEL or just on import side?
children = eval(mod + '_NAMES(''CHILDREN'');');
for i = 1:length(children)
    ichild = children{i};
    export_mod_dict(t, ichild, cond, "", mod, start);
    start = start + size_cond(ichild);
end
export_cond_dict(mod, cond, start);
if strlength(parent) == 0
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
if isempty(export_mod_out)
    export_mod_out = dictionary();
    export_mod_out('ALL') = VARS_Dict();
end
eval("export_mod_out('" + mod + "') = " + mod + "_Dict(cond, start);");
if strcmp(export_mod, mod) == 1
    global export_mod_file_base;
    if isempty(export_mod_file_base)
        export_mod_file_base = mod + "_";
    end
    if (t == 0)
        fname = export_mod_file_base + "init.txt";
    else
        fname = export_mod_file_base + "last.txt";
    end
    dump_vars(fname, export_mod_out);
    export_mod_enabled = 0;
end
end
