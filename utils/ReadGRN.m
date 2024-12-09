function out = ReadGRN(name)
order = [
"Glyma.19G046800", ...
"Glyma.08G165400", ...
"Glyma.08G165500", ...
"Glyma.04G015900", ...
"Glyma.03G185800", ...
"Glyma.10G268500", ...
"Glyma.07G142700", ...
"Glyma.10G293500", ...
"Glyma.11G226900", ...
"Glyma.15G136200", ...
"Glyma.05G025300", ...
"Glyma.01G026700", ...
"Glyma.19G089100", ...
"Glyma.06G094300", ...
"Glyma.10G210700", ...
"Glyma.04G011900", ...
"Glyma.17G015600", ...
"Glyma.19G022900", ...
"Glyma.07G028500", ...
"Glyma.16G168000", ...
"Glyma.19G017200", ...
"Glyma.08G044100", ...
"Glyma.06G323700", ...
"Glyma.10G086600", ...
"Glyma.11G169700", ...
"Glyma.09G015500", ...
"Glyma.15G012500", ...
"Glyma.08G097300", ...
"Glyma.06G017900", ...
"Glyma.08G302600", ...
"Glyma.09G255200", ...
"Glyma.04G051200", ...
"Glyma.13G222300", ...
"Glyma.13G204800", ...
"Glyma.10G042000", ...
"Glyma.01G095900", ...
"Glyma.09G024100", ...
];
InputF = importdata(name);
DataCol = InputF.textdata(1,:);
idx = find(ismember(DataCol, 'Ele:Amb' ));
names = InputF.textdata(:,1);
if (idx)
    idx = idx - 2;
    vals = InputF.data(:, idx);
else
    try
        map = ReadParam(name);
        names = keys(map);
        vals = values(map);
    catch ME
        InputF.textdata
        error("Could not determine the correct column to get expressions from");
    end
end
out = ReorderData(order, names, vals, 1.0);
end
