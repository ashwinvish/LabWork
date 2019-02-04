function [inputs] = inputsOfClass(class, df)

inputs = [];
for i= 1:size(class,2)
    [A,B] = SynapticPartners(class(i),1,df); % presynaptic partners
   % temp = A(A<1e5);
    inputs(1:size(A,1),i) = A;
end
end