function WriteTable(fname, names, values)
if strlength(fname) == 0
    return;
end
fileID = fopen(fname, 'w');
if length(names) > 0
    nameLine = strjoin(names, ",");
    fprintf(fileID, nameLine);
    fprintf(fileID,"\n");
end
valueStr = cell(size(values));
for i=1:length(values)
    valueStr{i} = num2str(values(i));
end
valueLine = strjoin(valueStr, ",");
fprintf(fileID, valueLine);
end
