% This function randomly draws 10% of current system and transfers it to a
% new fresh cell culture of 90% capacity.

function [new_cell, new_vir_out] = new_sample(cell, vir_out, sim_para)

global no_vir no_cu no_cv no_cd

%% Drawing 10% of current passage 
sample_vir = ceil(no_vir/10); % New sample size
sample_cell = ceil((no_cu + no_cv + no_cd)/10); % New sample size

vir_pos = randi(no_vir, sample_vir, 1); % Drawing random virus positions
cell_pos = randi((no_cu + no_cv + no_cd), sample_cell,1); % Drawing random cell positions


% Updating new virus profiles
new_vir_out.type = vir_out.type(vir_pos);
new_vir_out.profile = vir_out.profile(vir_pos,:);
no_vir = sample_vir;

%% New fresh cell culture of 90% capacity
no_new_cell = ceil(0.9*sim_para.i_no_cell);

% Updating new cell profiles
new_cell.type    = ones(no_new_cell,1);   % Initially all are uninfected cells
new_cell.nvir    = zeros(no_new_cell,1);  % Initially no virus
new_cell.vir_pop = repmat(struct('vir_in',struct('type',{0},'profile',{0})),no_new_cell,1);
new_cell.tau     = zeros(no_new_cell,1);  % Post infection time.(by normal virus only)

% Now adding old cells (10%)
for i=1:sample_cell
    new_cell.type(no_new_cell+i)    = cell.type(cell_pos(i));
    new_cell.nvir(no_new_cell+i)    = cell.nvir(cell_pos(i));
    new_cell.vir_pop(no_new_cell+i) = cell.vir_pop(cell_pos(i));
    new_cell.tau(no_new_cell+i)     = cell.tau(cell_pos(i));
end
    
return
