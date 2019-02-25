for i = 1:size(ContraRhombomere,2)
    A(i,:) = struct2array(ContraRhombomere(i));
end

index = 1
for i = 1:size(ContraRhombomere,2)
    if find(A(i,2:end)) == 1
        r3Contra(index) = A(i,1);
        index = index+1;
    end
end

index = 1
for i = 1:size(ContraRhombomere,2)
    if find(A(i,2:end)) == 2
        r4Contra(index) = A(i,1);
        index = index+1;
    end
end

index = 1
for i = 1:size(ContraRhombomere,2)
    if find(A(i,2:end)) == 3
        r5Contra(index) = A(i,1);
        index = index+1;
    end
end

index = 1
for i = 1:size(ContraRhombomere,2)
    if find(A(i,2:end)) == 4
        r6Contra(index) = A(i,1);
        index = index+1;
    end
end

index = 1
for i = 1:size(ContraRhombomere,2)
    if find(A(i,2:end)) == 5
        r7Contra(index) = A(i,1);
        index = index+1;
    end
end