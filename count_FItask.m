function [FI_f_bin,FI_f_sizes] = ...
    count_FItask(subG_f_bin,subG_f_sizes,att_inc_bound_f)
FIsubG=0;
FI_f_bin=zeros(size(subG_f_bin,2), 1);

for subGIdx=1:size(subG_f_sizes,2)
    LogIdx=find(subG_f_bin==subGIdx);
    f_test=LogIdx(1);
    if att_inc_bound_f(f_test)==true
        subG_f_bin(LogIdx)=0;
        subG_f_sizes(subGIdx)=0;
    else
        FIsubG=FIsubG+1;
        FI_f_bin(LogIdx)=FIsubG;
        FI_f_sizes(FIsubG)=subG_f_sizes(subGIdx);
    end
end