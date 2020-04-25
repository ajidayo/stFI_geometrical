function [LastVal] ...
    = execution_with_TMM_Explicit(TMM_Explicit,InitVal,number_of_steps)

%disp('checkpoint alpha')

variables_f_then_e=[InitVal.f; InitVal.e];

for k=1:number_of_steps
    variables_f_then_e=TMM_Explicit*variables_f_then_e;
end

LastVal = variables_f_then_e;

end

%disp('checkpoint bravo')
