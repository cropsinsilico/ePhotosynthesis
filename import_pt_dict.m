function out = import_pt_dict(mod, pt)
if (pt == "VEL")
    eval("names = " + mod + "_NAMES(pt);");
    eval("global " + mod + "_VEL;");
    eval("global " + mod + "_TIME_N;");
    eval("vels = " + mod + "_VEL;");
    eval("N = " + mod + "_TIME_N;");
    out = dictionary();
    if (N > 0)
        if (size(vels, 1) == N)
            out = dictionary(names, ...
                vels(N, 2:(2 + length(names) - 1)));
        elseif (size(vels, 2) == N)
            out = dictionary(names, ...
                vels(2:(2 + length(names) - 1), N)');
        end
    end
else
    eval("global " + mod + "_" + pt + "_DICT;");
    eval("out = " + mod + "_" + pt + "_DICT;");
end
end
