% This function distributes extracellular viruses into bins (We can
% consider different bin lengths also).
% Here we neglected viruses with length lesser than 10 kb, since tracking
% them in experiments are almost impossible.

function [bin, vir_range] = virus_range_new(vir_len)

global no_vir
vir_range = zeros(7,1); % Initialization
bin = [35 63 70 83 97 109 133.9];   % binpositions

vir_range(1) = length(find(vir_len>= 10 & vir_len <50));%/no_vir; %%$%#^%#^%^%
vir_range(2) = length(find(vir_len>= 50 & vir_len <65));%/no_vir;
vir_range(3) = length(find(vir_len>= 65 & vir_len <75));%/no_vir;
vir_range(4) = length(find(vir_len>= 75 & vir_len <85));%/no_vir;
vir_range(5) = length(find(vir_len>= 85 & vir_len <100));%/no_vir;
vir_range(6) = length(find(vir_len>= 100 & vir_len <125));%/no_vir;
vir_range(7) = length(find(vir_len>= 125 & vir_len <135));%/no_vir;

s =sum(vir_range);
for i=1:7
    vir_range(i) = vir_range(i)/ s;
end

return