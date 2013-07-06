function d = get_list_txt_files(dirname)
% LIST_TXT_FILES = GET_LIST_TXT_FILES(DIRNAME)
%
% reads only '*.txt' files from the DIRNAME, ignoring
% '_flt' and '_noflt' files 



% Keep current directory, change to dirname, read all txt files

% wd = cd;
% cd(dirname);
% fullfile removes the need to change directories [AL, 8.6.13]
d1 = dir(fullfile(dirname,'*.txt'));
d2 = dir(fullfile(dirname,'*_noflt.txt'));
d3 = dir(fullfile(dirname,'*_flt.txt'));

% the following loops remove _flt and _noflt.txt from the list
% of names. 

d = d1;
ind = false(1,length(d1));

for i = 1:length(d1)
    for j = 1:length(d2)
        if strcmp(d1(i).name,d2(j).name), ind(i) = true; end
    end
    for j = 1:length(d3)
        if strcmp(d1(i).name,d3(j).name), ind(i) = true; end
    end
end
d(ind) = [];