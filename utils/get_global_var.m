function value = get_global_var(name)
value = [];
try
    eval("global " + name + ";");
    eval("value = " + name + ";");
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:m_illegal_character')
        rethrow(ME);
    end
end
end
