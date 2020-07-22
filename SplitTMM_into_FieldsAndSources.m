function [TMM_Fields, TMM_Sources]   = SplitTMM_into_FieldsAndSources(TMM,Source,Num_of_Elem)
Sum=0;
for nth_source = 1:size(Source,2) 
    Sum=Sum+Source(nth_source).UpdNum;
end
STSourceNum = Sum;

TMM_Fields  = TMM(1:Num_of_Elem.SpP+Num_of_Elem.SpS,1:Num_of_Elem.SpP+Num_of_Elem.SpS);
TMM_Sources = TMM(1:Num_of_Elem.SpP+Num_of_Elem.SpS,...
    Num_of_Elem.SpP+Num_of_Elem.SpS+1:Num_of_Elem.SpP+Num_of_Elem.SpS+STSourceNum);


end