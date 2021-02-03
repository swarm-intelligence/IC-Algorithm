function [Out, nfe, fitness_result] = ICA_Sort(X,n_size,nfe,dim,max_nfe)
s = size(X,1);
F_result = zeros(1, s);
for i = 1:s
    F_result(1,i) = F3(X(i,:));
    nfe = nfe + 1;
    if(nfe >= max_nfe)
       break;
    end
end

[F_result, sorted_index] = sort(F_result,2);

sorted = zeros(n_size,dim);

for j=1:n_size
    sorted(j,:) = X(sorted_index(1,j),:);
end
Out = sorted;
fitness_result = F_result;


