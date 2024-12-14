function out = vardict2string(x)
out = "";
if isempty(x) || ~isConfigured(x)
    return;
end
modules = ["BF", "FIBF", "FI", "PR", "PS", "SUCS", "ALL"];
for imod = 1:length(modules)
    mod = modules(imod);
    if (~isKey(x, mod))
        continue
    end
    xmod = x(mod);
    if isempty(xmod) || ~isConfigured(xmod)
        continue;
    end
    if (mod == "ALL")
        param_types = keys(xmod);
    else
        param_types = ["MOD", "COND", "POOL", "RC"];
    end
    for ipt = 1:length(param_types)
        pt = param_types(ipt);
        if (isKey(xmod, pt))
            xpt = xmod(pt);
            if isempty(xpt) || ~isConfigured(xpt)
                continue;
            end
            vars = keys(xpt);
            for ivar = 1:length(vars)
                k = vars(ivar);
                xk = xpt(k);
                out = out + var2string(mod, pt, k, xk);
            end
        end
    end
end
end
