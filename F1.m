%Rastrigin
function [Out] = F1(X)
T1 = X.^2;
T2 = -10*cos(2*pi*X)+10;
Out = sum(T1+T2,2);
%D=[-5.12,5.12]