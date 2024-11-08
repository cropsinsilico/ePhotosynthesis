function values = get_global_values(names)
values = zeros(size(names));
for i = 1:length(names)
    iname = names{i};
    if iname(1:1) ~= "_"
        x = get_global_var(iname);
        if ~isempty(x)
            values(i) = x;
        end
    end
end
end
