function [list] = readintfile(name, nskip, cols)

% Deze functie leest tabgescheiden kolom text files
% de getallen worden naar float geconverteerd
% het gaat mis als de file andere dingen dan getallen bevat
% alleen Nan (not a number) en inf worden geaccepteerd

fid         = fopen(name,'rt');
for p=1:nskip
    fgetl(fid);
end
str         = fread(fid,inf,'*char');
st          = fclose(fid);
list        = sscanf(str','%f');

if nargin>=3
    assert(mod(length(list),cols)==0,'Number of columns in file not as expected or file not complete. Got %d elements from file',length(list))
    list = reshape(list,cols,length(list)/cols).';
end