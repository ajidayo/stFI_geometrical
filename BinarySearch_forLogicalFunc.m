function [Solution]=BinarySearch_forLogicalFunc(FuncHandle,FuncDomain,SearchEpsilon,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14)
WidthOfSolution=FuncDomain.to-FuncDomain.from;
% assuming that the function is true for x s.t. x<Solution, false for x s.t. x>Solution
Funcistrue=FuncHandle(WidthOfSolution/2+FuncDomain.from,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14);
FuncDomain.from = Funcistrue*(WidthOfSolution/2)+FuncDomain.from;
FuncDomain.to = not(Funcistrue)*(-WidthOfSolution/2)+FuncDomain.to;
disp(FuncDomain)

if WidthOfSolution<SearchEpsilon
    Solution=WidthOfSolution/2+FuncDomain.from;
else
    Solution=BinarySearch_forLogicalFunc(FuncHandle,FuncDomain,SearchEpsilon,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14);
end

end