%Schwefel 2.26
function [Out] = F3(X)
T1 = 418.9829.*size(X,2);
T2 = sum(X.*sin(sqrt(abs(X))),2);
Out = T1 - T2;
%D=[-500,500]