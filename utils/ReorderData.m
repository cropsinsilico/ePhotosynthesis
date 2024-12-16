function out = ReorderData(order, keys, values, default)
arguments
    order (:,1) string
    keys (:,1) string
    values (:,1) {mustBeNumeric}
    default (1,1) {mustBeNumeric} = 0
end
N = size(order, 1);
out = default * ones(N, 1);
for i = 1:N
    idRow = find(ismember(keys, order(i)));
    if idRow
        out(i, 1) = values(idRow, 1);
    end
end
end
