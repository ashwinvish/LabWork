temp = [];
intPartners = [PrePartners, PostPartners];
for i = 1: size(intPartners,2)
    temp = [temp ; intPartners{i}];
end
intPartners = unique(temp(temp<1e5)); 
IntConnMatrixPre = zeros(1,size(intPartners,1));
for i = 1:1%size(intPartners,1) % pre
        tempPrePartner =  df.presyn_seg(df.postsyn_seg==77803);
        tempPrePartner = tempPrePartner(tempPrePartner<1e5);
        if ~tempPrePartner == 0
            %i
            [N,edges] = histc(tempPrePartner, unique(tempPrePartner));
            [a,b] = intersect(intPartners,unique(tempPrePartner));
            for j = 1:size(b,1)
             IntConnMatrixPre(i,b(j)) = N(find(a(j)==unique(tempPrePartner)));
            end
        else 
            continue;
        end
        IntConnMatrixPre(i,i) = 0;
        clear tempPrePartner;      
end
