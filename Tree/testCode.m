tic
for i = 1:imHeight
    for j = 1:imNumber-1
        SideImageXZ(j,:,i) = FinalImage(i,:,j);
    end
end
toc
