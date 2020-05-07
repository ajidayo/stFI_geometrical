function [s_or_t,taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task]...
    = semi_add_dep_edge(p,taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task,p_is_init)
if taskIdx_for_p(p)==0
    taskIdx_newest=taskIdx_newest+1;
    s_or_t=taskIdx_newest;
    taskIdx_for_p(p)=taskIdx_newest;
    task(taskIdx_newest).p_tgt=p;
    if p_is_init==true
        task(taskIdx_newest).typ="InitVal";
    else 
        task(taskIdx_newest).typ="stFI";
    end
else
    s_or_t=taskIdx_for_p(p);
end
end