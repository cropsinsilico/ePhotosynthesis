function out = size_cond(mod)
    out = length(eval(mod + "_NAMES('COND');"));
    children = eval(mod + "_NAMES('CHILDREN');");
    for i = 1:length(children)
        ichild = children{i};
        out = out + size_cond(ichild);
    end
end
