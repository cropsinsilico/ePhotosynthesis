function fname = dump_vars(fname, vars)
fileID = fopen(fname,'w');


varStr = vardict2string(vars);
fprintf(fileID, '%s', varStr);

end