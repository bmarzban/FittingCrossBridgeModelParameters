function [choosed] =  roulettewheel2(relative_matrix,alimit,aa)

rand_prob = ceil(rand*1000);
for i=1:aa

if floor(alimit(i)/rand_prob)>=1
    choosed = i;
break

end
  
end



