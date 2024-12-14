function out = compare_floats(actual, expected, reltol, abstol, fmt)
    arguments
        actual (:, :) double
        expected (:, :) double
        reltol (1, 1) double = 1.0e-03
        abstol (1, 1) double = 1.0e-03
        fmt (1, 1) string = ""
    end
    if strlength(fmt) == 0
        nprec = floor(max(log10(1/abstol), log10(1/reltol))) + 1;
        width = nprec + 6;
        fmt = ['% ' num2str(width) '.' num2str(nprec) 'f '];
    end
    reldiff = abs((actual - expected) ./ actual);
    absdiff = abs(actual - expected);
    out = 1;
    if ~all(reldiff <= reltol)
        out = 0;
        fprintf("Relative difference (%s) exceeds tolerance (%s)\n", ...
                sprintf(fmt, reldiff), sprintf(fmt, reltol));
    end
    if ~all(absdiff <= abstol)
        out = 0;
        fprintf("Absolute difference (%s) exceeds tolerance (%s)\n", ...
                sprintf(fmt, absdiff), sprintf(fmt, abstol));
    end
    if (out == 0)
        fprintf("Actual      %s\nExpected    %s\n", ...
                sprintf(fmt, actual), sprintf(fmt, expected));
    end
end
