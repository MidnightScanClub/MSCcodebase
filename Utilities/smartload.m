function out = smartload(matfile)

out = load(matfile);
names = fieldnames(out);
out = eval(['out.' names{1}]);