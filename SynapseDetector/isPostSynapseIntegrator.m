function [logic] = isPostSynapseIntegrator(cellID,df)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

confirmedALX = [76181 76201 76187 76184 76192 76197 ];
confirmedDBX = [76199 76200 76182 76183 76185 76186 76189 76191 76188 ];
confirmedBARHL = [76198 76190 76193 76194 76195 76196 ];
allIntegrators  = [confirmedALX,confirmedDBX,confirmedBARHL];

logic = zeros(size(cellID));
for i = 1:numel(cellID)
    [A,~] = SynapticPartners(cellID(i),2,df);
    valid = ismember(A,allIntegrators);
    if sum(valid)>0
        logic(i) = 1;
    end
    clear A;
    clear valid;
end


end

