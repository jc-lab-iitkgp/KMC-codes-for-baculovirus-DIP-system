% This function selects the cell and virus with uniform randomness.
% This program executes after the event is identified.
function [cellwhich, viruswhich] = cell_vir_select(cell, eventwhich, rate_per_cv, sys_para)

global no_vir 

% Identification of locations in the cell structure where
% uninfected, DIP infected and virus infected cells are present. 
CU_pos = find(cell.type ==1);  %Positions of CU
CD_pos = find(cell.type ==2);  %Positions of CD
CV_pos = find(cell.type ==3);  %Positions of CV

switch eventwhich
    
    case 1  % Outside virus entering in CU
        cellwhich  = datasample(CU_pos,1);
        viruswhich = randi(no_vir);
        
    case 2  % Outside virus infecitng into CV
        cellwhich  = datasample(CV_pos,1);
        viruswhich = randi(no_vir);
        
    case 3  % Outside virus infecting into CD
        cellwhich  = datasample(CD_pos,1);
        viruswhich = randi(no_vir);
 
    case 4  % Replication of uninfected cells
        cellwhich  = datasample(CU_pos,1);
        viruswhich = [];
       
    case 5  % Replication of DIP infected cells
        cellwhich  = datasample(CD_pos,1);
        viruswhich = [];
        
    case 6  % Cell death
        cum_sum_death_rate = cumsum(rate_per_cv.death);
        random_no          = rand* cum_sum_death_rate(end); 
        ind                = find(cum_sum_death_rate >= random_no,1);  
        cellwhich          = CV_pos(ind);
        viruswhich         = [];
        
    case 7  % virus replication
        cum_sum_vrepli_rate = cumsum(rate_per_cv.vrepli);        
        random_no = rand* cum_sum_vrepli_rate(end); 
        ind = find(cum_sum_vrepli_rate >= random_no,1); 
        cellwhich  = CV_pos(ind); % Selected cell
        
        % %% THIS PART TO USE IF ALL VIRUSES HAVE EQUAL PROB. TO GET SELECTED FOR VIRUS REPLICATION
        viruswhich = randi(cell.nvir(cellwhich)); % Selected virus
        
        % %% THIS PART TO USE IF SMALLER VIRUSES HAVE MORE PROB. TO GET SELECTED FOR VIRUS REPLICATION
%         for i=1: cell.nvir(cellwhich)
%             inv_len(i) = 1 /(cell.vir_pop(cellwhich).vir_in(i).profile * sys_para.hr_kb') ; % Inverse of lengths of viruses inside chosen cell
%         end
%         cum_sum_inv_len = cumsum(inv_len);
%         viruswhich = find(cum_sum_inv_len >= rand* cum_sum_inv_len(end), 1); % Selected virus
                
    case 8  % Virus release 
        cum_sum_vrel_rate = cumsum(rate_per_cv.vrel);
        random_no = rand* cum_sum_vrel_rate(end); 
        ind = find(cum_sum_vrel_rate >= random_no,1);
        cellwhich  = CV_pos(ind); % Selected cell
        viruswhich = randi(cell.nvir(cellwhich)); % Selected virus     
end   
        
return