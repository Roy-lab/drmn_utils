function quantile_normalize(ifname, outfname, dolog, groups)
%% Does log transform and quantile normalization on an input gene exp file
% If plotting filename prefix provided, will make three output plots: 
% two containing quantile values and one showing gene expression values
% before/after.
%
% Args:
%   ifname: Gene expression file with "Gene" in top left; nodata values
%           represented with NaN ONLY
%   outfname: name for output expression file
%   plotpref (optional): name for output plot files
% Original by SR; Updated by DC with plotting

% check how many args were provided
% if no plot filename prefix, then don't make plots
switch nargin
    case 2
        fprintf(2,'Bad argument list\n');
        return
    case 3
        plotpref='';
end
        
d=importdata(ifname,'\t')
%order=[22:24 1:3 13:15 4:6 10:12 16:21 7:9]
data=d.data(:,:);
if (strcmp(dolog,'doLog'))
    ldata=log(data+1);
else
    ldata=data;
end

normdata=quantilenorm(ldata);
size(normdata)
% if plot filename given, will make plots of before/after

% print out normalized log data
fid=fopen(outfname,'w');
for r=1:size(normdata,1)
    fprintf(fid,'%s',d.textdata{r});
    for c=1:size(normdata,2)
        fprintf(fid,'\t%f',normdata(r,c));
    end
    fprintf(fid,'\n');
end
fclose(fid);
end

