function out = ReadParam(name)
data = readtable(name);
keys = string(data{:, 1});
values = data{:, 2};
out = dictionary(keys, values);
end
