function [ popNew ] = crossover( popPrev, crossProb )
%CROSSOVER Generate the new population by crossing over the paramters
% output:
% popNew: the new generated population
% input:
% popPrev: the previous population
% crossProb: the probability of crossing over
[new, sortIndex] = sort(rand(size(popPrev ,1) ,1));
popSorted = popPrev(sortIndex, :);
pairsNum = size(popSorted, 1)/2;
chromosomeEachSizes = size(popSorted, 2);
parisToCross = rand(pairsNum, 1) < crossProb;
pointsToCross = parisToCross.*randi([1,chromosomeEachSizes],[pairsNum, 1]);
Cross_over_rnd = rand(pairsNum,1);

for i=1:pairsNum
    if pointsToCross(i)==0
      
        popNew([2*i-1,2*i],:)=[popSorted([2*i-1,2*i],1:pointsToCross(i)), popSorted([2*i,2*i-1],pointsToCross(i)+1:chromosomeEachSizes)];
    else
%         popNew([2*i-1,2*i],:)=[popSorted([2*i-1,2*i],1:pointsToCross(i))*Cross_over_rnd(i), popSorted([2*i,2*i-1],pointsToCross(i)+1:chromosomeEachSizes)*(1-Cross_over_rnd(i))];
        
        popNew([2*i-1,2*i],:)=[popSorted(2*i-1,1:pointsToCross(i))*Cross_over_rnd(i), popSorted(2*i,pointsToCross(i)+1:chromosomeEachSizes);
                               popSorted(2*i,1:pointsToCross(i))*(1-Cross_over_rnd(i)), popSorted(2*i-1,pointsToCross(i)+1:chromosomeEachSizes)];
       

    end
    
end
end