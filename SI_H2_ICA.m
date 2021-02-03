clc;
clear;
close all;

%Number of Population & Dimensions
pop_size = 100;
dim = 2;

%Iteration Condition
max_iteration = 220;

%Domain of Benchmarks
from = -500;
to = -1*from;

%Results for n Times Execution
num_of_result = 5;
%Columns of total Result : dim,(gbest_fitness),(time)
total_result = zeros(num_of_result,dim+2);

for n=1:num_of_result
    
    tic;
    
    %Nfe Condition
    nfe = 0;
    max_nfe = 20000;
    
    %Initializing Random Population
    Countries = unifrnd(from,to,[pop_size dim]);
    
    %Number of the Empires
    n_empire = ceil(0.1 * pop_size);
    Empires = zeros(n_empire,dim+1);
    
    %Number of Colonies
    n_col = pop_size - n_empire;
    Colonies = zeros(n_col,dim+2);
    
    %Minimum Colonies for each Empire
    min_col = 3;
    
    %Revolution_Rate
    rev_rate = 0.05;
    
    %Calculating Fitness & Sorting
    [sorted_Countries, nfe, F_result] = ICA_Sort(Countries,pop_size,nfe,dim,max_nfe);
    
    %Initializing Empires
    Empires(:,1:dim) = sorted_Countries(1:n_empire,:);
    Empires(:,dim+1) = F_result(1,1:n_empire);
    
    %Initializing Colonies
    Colonies(:,1:dim) = sorted_Countries(n_empire+1:end,:);
    Colonies(:,dim+1) = F_result(1,n_empire+1:end);
    
    %Normalized Fitness of each Empire
    emptyE = zeros(1,n_empire);
    cn = emptyE;
    Cn = emptyE;
    cn(1,:) = Empires(:,dim+1);
    Cn(1,:) = max(cn) - cn(1,:);
    
    %Relative Power of each Empire
    pn = emptyE;
    pn(1,:) = abs(Cn(1,:) / sum(Cn));
    
    %Number of Primitive Colonies for each Empire
    NCn = emptyE;
    NCn(1,:) = min_col;
    NCn(1,:) = NCn(1,:) + (round( pn(1,:) * (n_col - (n_empire*min_col))) );
    
    %Shuffling the Colonies
    randx = randperm(n_col);
    Colonies = Colonies(randx,:);
    
    %Assign Colonies to the Empires
    j = 1;
    k = NCn(1,1);
    for i = 1:n_empire
        Colonies(j:k,dim+2) = i;
        j = j + NCn(1,i);
        if i+1 > n_empire
            break;
        end
        k = k + NCn(1,i+1);
    end
    
    %Control the Assignation Term
    if size(Colonies,1) > n_col
        Colonies = Colonies(1:n_col,:);
    end
    %Control the Assignation Term
    for i = 1:n_col
        if Colonies(i,dim+2) == 0
            rand_empire = randi(n_empire);
            Colonies(i,dim+2) = rand_empire;
        end
    end
    
    %Starting the Main Loop
    for iteration = 1:max_iteration
        
        %Compute Distance Between the Empires and their Colonies (Euclidean)
        distance = zeros(n_col,dim);
        for i = 1:n_empire
            for j = 1:n_col
                if Colonies(j,dim+2) == i
                    distance(j,:) = sqrt(sum((Empires(i,1:dim) - Colonies(j,1:dim)).^2));
                end
            end
        end
        
        %Move the Colonies
        beta = 2;
        teta = unifrnd(-pi/4,pi/4,[n_col dim]);
        x = zeros(n_col,dim);
        x(:,:) = unifrnd(0,beta*distance(:,:));
        X = teta .* x;
        Colonies(:,1:dim) = X + Colonies(:,1:dim);
        
        %Control the Domain of the Colonies
        for i = 1:n_col
            for j = 1:dim
                r_c = rand;
                if (Colonies(i,j) > to)
                    Colonies(i,j) = to - r_c;
                elseif (Colonies(i,j) < from)
                    Colonies(i,j) = from + r_c;
                end
            end
        end
        
        %Fitness of New Colonies
        Colonies(:,dim+1) = F3(Colonies(:,1:dim));
        nfe = nfe + n_col;
		
		%Termination Condition
        if(nfe >= max_nfe)
            break;
        end
        
        %Exchange the Positions of Empire and Colony
        %(if Fitness of Colony is better than Empire)
        for i = 1:n_empire
            for j = 1:n_col
                if Colonies(j,dim+2) == i
                    if Colonies(j,dim+1) < Empires(i,dim+1)
                        temp = Empires(i,1:dim+1);
                        Empires(i,1:dim+1) = Colonies(j,1:dim+1);
                        Colonies(j,1:dim+1) = temp;
                    end
                end
            end
            %Calculate the Average Fitness of the Colonies of each Empire
            ind = find(Colonies(:,dim+2)==i);
            if(isempty(ind) == 0)
                Empires(i,dim+2) = mean(Colonies(ind(1):ind(end),dim+1));
            end
        end
        
        %The Whole Power of each Empire
        zeta = 0.05;
        TCn = Empires(:,dim+1) + (zeta * Empires(:,dim+2));
        
        %Normalized The Whole Power of each Empire
        NTCn = max(TCn) - TCn;
        
        %Probability of Takeover for each Empire
        Ppn = abs(NTCn / sum(NTCn));
        
        %Colonial Competiton
        if (size(Empires,1) > 1)
            R = rand(n_empire,1);
            D = Ppn - R;
            weak_emp_ind = find(NTCn==min(NTCn)); %Minimization
            weak_colonies = sortrows(Colonies(Colonies(:,dim+2)==weak_emp_ind,:),-(dim+1));
            colonies_takeover_size = ceil(size(weak_colonies,1) * 0.5);
            takeover_ind = find(D==max(D));
            for i = 1:colonies_takeover_size
                weakest_colonies_ind = find((Colonies(:,dim+2)==weak_emp_ind) & (Colonies(:,dim+1) == weak_colonies(i,dim+1)),1);
                Colonies(weakest_colonies_ind,dim+2) = takeover_ind;
            end
        end
        
        %Fall of the Weakest Empire
        if (size(weak_colonies,1) - colonies_takeover_size) == 0
            if n_empire - 1 ~= 0
                n_col = n_col + 1;
                n_empire = n_empire - 1;
                Colonies(n_col,1:dim+1) = Empires(weak_emp_ind,1:dim+1);
                rand_empire = randi(n_empire);
                Colonies(n_col,dim+2) = rand_empire;
                Empires = removerows(Empires,'ind',weak_emp_ind);
                %Reindex
                for index = n_empire+1:-1:weak_emp_ind
                    if index - 1 ~= 0
                        Colonies(Colonies(:,dim+2)==index,dim+2) = index-1;
                    end
                end
            end
        end
        
        %Revolution Process
        y = unifrnd(from,to,[floor(rev_rate*n_col) dim]);
        for i = 1:size(y,1)
                rand_rev = randi(n_col);
                Colonies(rand_rev,1:dim) = Colonies(rand_rev,1:dim) + y(i,:);
        end
        
    end
    
    total_result(n,1) = toc;
    %total_result(n,2) = min(Empires(:,dim+1));
    %t_r_d = find(min(Empires(:,dim+1)));
    %total_result(n,3:end) = Empires(t_r_d,1:dim);
    total_result(n,2) = Empires(1,dim+1);
    total_result(n,3:end) = Empires(1,1:dim);
end

min_fitness = min(total_result(:,2));
max_fitness = max(total_result(:,2));
mean_fitness = mean(total_result(:,2));
std_fitness = std(total_result(:,2));
mean_time = mean(total_result(:,1));

disp(strcat('Popsize:', num2str(pop_size), ', Dim:', num2str(dim)));
disp(strcat('mean fitness: ', num2str(mean_fitness)));
disp(strcat('max fitness: ', num2str(max_fitness)));
disp(strcat('min fitness: ', num2str(min_fitness)));
disp(strcat('std fitness: ', num2str(std_fitness)));
disp(strcat('mean time: ', num2str(mean_time)));