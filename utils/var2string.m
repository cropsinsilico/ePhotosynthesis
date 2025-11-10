function out = var2string(mod, pt, name, var)
prefix = mod + "::" + pt + "::" + name;
fmt = "%-35s  %g\n";
out = sprintf(fmt, prefix, var);
end