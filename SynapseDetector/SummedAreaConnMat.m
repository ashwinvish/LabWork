% connmatrix with summed area

load AllCells;
LoadDataFrame;

connMatSummed = zeros(size(AllCells));
for i = 1:size(AllCells,1)
    temp = SynapticPartners(AllCells(i),1,df); % all presynaptic partners
    [N]= histc(temp,AllCells);
    index = find(N) ; % all positive entires i.e all partners that make synapses
    for j = 1:length(index)
    sz = df.size(df.presyn_segid == AllCells(index(j)) & df.postsyn_segid == AllCells(i));
    connMatSummed(i,index(j)) = sum(sz);
    clear sz;
    end  
    clear temp;
    connMatSummed(i,i) = 0;
end
ConnMatSummed = connMatSummed;

save('ConnMatSummed.mat','ConnMatSummed');