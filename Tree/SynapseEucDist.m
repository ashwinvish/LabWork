
index = 1;
for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        for jj = 1:length(allPreSynapse{i})
            AlxPreEuc(index,jj) = pdist([CellSoma(i,:); allPreSynapse{i}(jj,:)]);
        end
        for kk = 1:length(allPost{i})
            AlxPostEuc(index,kk) = pdist([CellSoma(i,:); allPost{i}(kk,:)]);
        end
        index = index +1;
    end
    
end

index = 1;
 for i = 1:numel(cellIDs)
     if ismember(cellIDs{i}, cellIDsTrans) ==1
        for jj = 1:length(allPreSynapse{i})
            TransPreEuc(index,jj) = pdist([CellSoma(i,:); allPreSynapse{i}(jj,:)]);
        end
        for kk = 1:length(allPost{i})
            TransPostEuc(index,kk) = pdist([CellSoma(i,:); allPost{i}(kk,:)]);
        end
        index = index +1;
     end  
 end

 index = 1;
 for i = 1:numel(cellIDs)
     if ismember(cellIDs{i}, cellIDsDbx) ==1
        for jj = 1:length(allPreSynapse{i})
            DbxPreEuc(index,jj) = pdist([CellSoma(i,:); allPreSynapse{i}(jj,:)]);
        end
        for kk = 1:length(allPost{i})
            DbxPostEuc(index,kk) = pdist([CellSoma(i,:); allPost{i}(kk,:)]);
        end
        index = index +1;
     end  
 end

 index = 1;
 for i = 1:numel(cellIDs)
     if ismember(cellIDs{i}, cellIDsL) ==1
        for jj = 1:length(allPreSynapse{i})
            BarhlPreEuc(index,jj) = pdist([CellSoma(i,:); allPreSynapse{i}(jj,:)]);
        end
        for kk = 1:length(allPost{i})
            BarhlPostEuc(index,kk) = pdist([CellSoma(i,:); allPost{i}(kk,:)]);
        end
        index = index +1;
     end  
 end

 
 subplot(4,1,1)
 histogram(nonzeros(AlxPreEuc)/1000,'BinWidth',2,'FaceColor',[0,0.8,0]);
 hold on;
 histogram(nonzeros(AlxPostEuc)/1000,'BinWidth',2,'FaceColor',[0.9,0,0]);

 
 subplot(4,1,2);
 histogram(nonzeros(TransPreEuc)/1000,'BinWidth',2,'FaceColor',[0,0.8,0]);
 hold on;
 histogram(nonzeros(TransPostEuc)/1000,'BinWidth',2,'FaceColor',[0.9,0,0]);

 
 subplot(4,1,3);
 histogram(nonzeros(DbxPostEuc)/1000,'BinWidth',2,'FaceColor',[0.9,0,0]);

 
 subplot(4,1,4);
 histogram(nonzeros(BarhlPostEuc)/1000,'BinWidth',2,'FaceColor',[0.9,0,0]);

 
 
 