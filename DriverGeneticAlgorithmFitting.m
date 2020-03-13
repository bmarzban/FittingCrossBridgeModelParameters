% Optimization Driver/ Sensitivity analysis
clear; close all; clc;
% define the variables to be tuned and their corresponding numbers as
% following: 
nvar = 15; 
Mutation_rate = 0.25;
Cross_over_rate = 0.7;
Population_num = 40; %%%even numbers
Max_gen = 100;
Elitism_rate = 0.2;
% plotting while the optimization runs
flag_plot = 0;
%%%%%%%Producing the first Generetion%%%%%%%%
Initial_parent = rand(Population_num,nvar);
% Initial_parent = ones(Population_num,nvar);
%% Adding upper band and lower band to the parameters value
variables_min_rang = ones(1,nvar) - 0.5; %[min_varizbale_1 min_variable_2 ...]
variables_max_rang = ones(1,nvar) + 1;   %[max_varizbale_1 max_variable_2 ...]
initial_population = zeros(Population_num,nvar);
% % 
% variables_min_rang(10:11) = [0.5];
variables_max_rang(10:11) =[7];
variables_min_rang(12:14) = [0.99];
variables_max_rang(12:14) =[1.01];
variables_min_rang(4:8) = [0.2];
variables_max_rang(4:8) =[10];
% good_individual_3 = [1.8685    0.9927    1.9746*1    1.0887*0.99    1.7791    1.5887    1.6376    0.5657    1.8738    1.9557    1.6768    1.7204];%139
% good_individual_1 = [1.18783129356751*1.43,0.755063156386634,0.608302745025594,2.28505848337249,2.83554452432431,2.36363384319356,0.908906573915146,1.81657209630941,1.24834053130923,2.0,1.14406033829670,1.00373905093119,0.938658804993310,1.60483521822591 ,2 ,2 ]
% good_individual_1 = [ 0.9774*1.32    1.5817    0.9128*1.4    1.6309    1.4893    0.6127    1.5998*0.9    1.4116*1.12    0.5402    1.8157    1.7444    1.0733    1.1222    1.7827    2 0.6905]


for i=1:nvar
    
      a = variables_min_rang (i);
      b = variables_max_rang (i);
      initial_population ( : , i) = a + ( b - a ) * Initial_parent( : , i );

end
% load previos simulations
% load run2.mat parentTemp2
% initial_population(1:2,:)= parentTemp2(end-1:end,:);
% % % 
%  initial_population (6,:) = good_individual_1;
% % initial_population (7,:) = good_individual_2;
% initial_population (8,:) = good_individual_3;

Data.freq =[2 3 4 5 6 7 8 9 10 4 4 4];
% Set temperature fot the experiment environment
Data.TmpC = 37.5; % centigrade 
tic
%Initialization the function value
f = zeros(Population_num,1);
f_min = zeros(Population_num,1);
f_max = zeros(Population_num,1);
%Using the parfor caculate objective function 
for z = 1:Population_num
    tic
 try 
    [f(z)] = objective_fun_rat_data(initial_population(z,:) , Data, flag_plot);
 catch
                disp '** Error Catch11'
 end
    toc1 = toc;
end   
% sort objective function values and also sort the variables (initial population) based on the corresponding function values 

[f,ind] = sort(f);
ParentSorted = initial_population(ind,:);
Elitism_num = floor(Elitism_rate*Population_num);
Mating_pool_num = ceil(Population_num-Elitism_num);
Maxf = -inf;
f_min(1) = min(f);
f_max(1) = max(f); %best f
selected = zeros(1,Mating_pool_num);
parentTemp2 = zeros(Population_num,nvar);
s = 0;
ParentTemp1 = ParentSorted;
Max_each_iteration = zeros(Max_gen,1);
Gen = 0;

while Gen < Max_gen
%    try
  %%%%caculating the relative fitness%%%%%%%%
  % Shift the objective functions values to have a positive value
  Modified_fitness = abs(min(f))+abs(max(f))+f;
  F_sum = sum(Modified_fitness);
  Relative_fitness = (Modified_fitness/F_sum);
  %%%%%Use fitness-proportionate selection (Roulette Wheel) to form the Mating Pool%%%%%
  alimit = ceil(Relative_fitness*1000);
  alimit = cumsum(alimit);
  for i = 1:Mating_pool_num
        selected(i) = roulettewheel2(Relative_fitness,alimit,Population_num);
  end
  
for n = 1 : Mating_pool_num
    ParentTemp1(n,:) = ParentSorted(selected(n),:);
end
parentTemp2 = ParentTemp1;
parentTemp2(1:Mating_pool_num,:) = crossover(ParentTemp1(1:Mating_pool_num,:),Cross_over_rate);
Mutation_rate_rnd = rand(Population_num,nvar);
Mutation_rnd = rand(Population_num,nvar);
% make sure that the parameters are actually in the range of upper bound
% and lower band after the cross over

   for i = 1 : Population_num-Elitism_num
       for j = 1 : nvar
           if Mutation_rate_rnd(i,j) < Mutation_rate
               a = variables_min_rang(j);
               b = variables_max_rang(j);
               parentTemp2(i,j)= a+(b-a) * Mutation_rnd(i,j);
           end
       end
   end
     

for h = 1 : Population_num
    % make sure that the parameters are actually in the range of upper bound
    % and lower band after the cross over
    temp_min_index = (parentTemp2 ( h , : )  >  variables_min_rang);
    parentTemp2 ( h , : ) = parentTemp2 ( h , : ).*temp_min_index + variables_min_rang .*(1- temp_min_index);
    temp_max_index = (parentTemp2 ( h , : )  <  variables_max_rang);
    parentTemp2 ( h , : ) = parentTemp2 ( h , : ).*temp_max_index + variables_max_rang .*(1- temp_max_index);
    
    [f(h)] = objective_fun_rat_data (parentTemp2 ( h , : ) , Data, flag_plot);
end   

  [ f , ind ] = sort (f);

 parentTemp2 = parentTemp2 (ind , :);

 if max(f) > Maxf
     s = s +1
 Maxf(s) = max(f)
 Max_X = parentTemp2(Population_num,:)
 end
 
 if max(f) > 0.3
     Gen = Max_gen;
      break
  end
ParentSorted = parentTemp2;
ParentTemp1 = ParentSorted;
Gen = Gen + 1
Max_each_iteration(Gen) = max(f);
%   catch
%             disp '** Error Catch'
%    end
  save run2.mat
end

iteration = [ 1 : 1 : Gen];
disp ( Maxf(s))
disp ( Max_X)

figure(9)
clf
plot(iteration,Max_each_iteration,'k','LineWidth',1.5)
grid on;
xlabel('Iteration','FontSize',14,'FontName','Times New Roman');
ylabel('Cost Function','FontSize',14,'FontName','Times New Roman');
legend('GA')
save run2.mat
toc