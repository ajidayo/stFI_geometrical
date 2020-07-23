function DoFs_FacesThenEdges = TimeMarch(Num_of_Steps,Time,cdt,TMM_Fields,TMM_Sources,DoFs_FacesThenEdges,Source)
Sum = 0;
for SourceIdx = 1:size(Source,2)
    Sum = Sum + Source(SourceIdx).UpdNum;
end
DoF_Source = zeros(Sum,1); 
for StepCount = 1:Num_of_Steps
    %StepCount
    for SourceIdx = 1:size(Source,2)
        UpdNum_Source       = Source(SourceIdx).UpdNum;
        WaveformFunction    = Source(SourceIdx).WaveformFunctionHandle;
        Area_TargetDualFace = Source(SourceIdx).Area_TargetDualFace;
        WaveformSign = Source(SourceIdx).WaveformSign;
        for CurrentTimeSection = 1:UpdNum_Source
            Working_at_Time = Time + (CurrentTimeSection-1)*cdt/UpdNum_Source;
            DoF_Source(Source(SourceIdx).FirstST_SourceIdx-1+CurrentTimeSection) ...
                = WaveformSign*WaveformFunction(Working_at_Time)*Area_TargetDualFace;
        end
    end
    DoFs_FacesThenEdges = TMM_Fields*DoFs_FacesThenEdges +TMM_Sources*DoF_Source;
    Time = Time + cdt;
end
end
