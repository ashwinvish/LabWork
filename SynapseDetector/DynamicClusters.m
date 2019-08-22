
fname = 'neuron_groups.json'; 
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

for i = 1:size(val,1)
    valid(i).Saccade= isSaccade(val{i});
    valid(i).Vestibular= isVestibular(val{i});
    valid(i).Integrator= isIntegrator(val{i});
    valid(i).Contra= isContra(val{i});
    valid(i).Motor = isAbducens(val{i});
end
    
for i = 1:size(val,1)
    subplot(4,7,i)
    imagesc(valid(i).Motor);
    %daspect([1,1,1]);
    title(val(i,2));
end

