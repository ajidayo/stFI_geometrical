F_Bz_symmetrical_error=zeros(MeshParam.Size_X,MeshParam.Size_Y);

for j=1:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        Temp=F(i,j).Bz-F(j,i).Bz.';
        UnSymFlag=any(any(find(Temp>EPSILON)));
        if UnSymFlag
            disp(['unsymmetric found: i=',num2str(i),', j=',num2str(j)])
        end
        F_Bz_symmetrical_error(i,j)=UnSymFlag;
    end
end
UnSymFlag=any(any(find(F_Bz_symmetrical_error>EPSILON)));
if UnSymFlag
    disp(['unsymmetric found: i=',num2str(i),', j=',num2str(j)])
end
