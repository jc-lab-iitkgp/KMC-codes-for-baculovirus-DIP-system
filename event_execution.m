% This function executes the chosen event and update the system (cell and
% virus profiles).

function [cell, vir_out] = event_execution(eventwhich, cellwhich, viruswhich, cell, vir_out, delta_t, sys_para)

global no_cu no_cd no_cv no_vir no_dip;  

switch eventwhich
   
    case 1  % Infection of uninfected cells
        no_cu = no_cu -1;   % No of uninfected cells reduced by 1
        no_vir = no_vir -1; % No of extracelluar vir reduced by 1 
        
        cell.nvir(cellwhich) = 1; % No of vir inside cell increases from 0 to 1
        %Copying vir profile into cell
        cell.vir_pop(cellwhich).vir_in(1).type = vir_out.type(viruswhich);
        cell.vir_pop(cellwhich).vir_in(1).profile = vir_out.profile(viruswhich,:); 
        
        if vir_out.type(viruswhich) == 2 % DIP infection
            no_dip = no_dip -1; % No of DIP outside decreases by 1
            no_cd = no_cd +1;   % No of CD increases by 1
            cell.type(cellwhich) = 2; % This cell becomes DIP infected (CD)    
        else   % normal virus(NV) infection
            no_cv = no_cv +1;   % No of CV increases by 1
            cell.type(cellwhich) = 3; % This cell becomes CV
        end
        
        %Outside virus delete
        vir_out.type(viruswhich) =[];
        vir_out.profile(viruswhich,:) =[];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    case 2  % Infection of virus infected cells
        no_vir = no_vir -1; % No of extracelluar vir reduced by 1 
        
        cell.nvir(cellwhich) = cell.nvir(cellwhich) +1;  % No of vir inside cell increases
        %Copying vir profile into cell
        cell.vir_pop(cellwhich).vir_in(cell.nvir(cellwhich)).type = vir_out.type(viruswhich);
        cell.vir_pop(cellwhich).vir_in(cell.nvir(cellwhich)).profile = vir_out.profile(viruswhich,:); 
        
        if vir_out.type(viruswhich) == 2 % DIP infection
            no_dip = no_dip -1;  % No of DIP outside decreases by 1
        end
        
        %Outside virus delete
        vir_out.profile(viruswhich,:) =[];
        vir_out.type(viruswhich) =[];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    case 3  % Infection of DIP infected cells
        no_vir = no_vir -1; % No of extracelluar vir reduced by 1
        
        cell.nvir(cellwhich) = cell.nvir(cellwhich) +1;  % No of vir inside cell increases
        %Copying vir profile into cell
        cell.vir_pop(cellwhich).vir_in(cell.nvir(cellwhich)).type = vir_out.type(viruswhich);
        cell.vir_pop(cellwhich).vir_in(cell.nvir(cellwhich)).profile = vir_out.profile(viruswhich,:); 
        
        if vir_out.type(viruswhich) == 2 % DIP infection
            no_dip = no_dip -1;  % No of DIP outside decreases by 1
        else   % normal virus(NV) infection
            no_cv = no_cv +1;  % No of CV increases by 1
            no_cd = no_cd -1;  % No of CD decreases by 1
            cell.type(cellwhich) = 3;
        end
        
        %Outside virus delete
        vir_out.profile(viruswhich,:) =[];
        vir_out.type(viruswhich) =[];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    case 4  % Replication of uninfected cells
        no_cu   = no_cu +1;      % No of CU increases by 1
        no_cell = no_cu + no_cv + no_cd;    % Total No of cells 
        % cell profile of new uninfected cell
        cell.type(no_cell)    = 1; 
        cell.nvir(no_cell)    = 0;  
        cell.vir_pop(no_cell) = struct('vir_in',struct('type',{0},'profile',{0}));
        cell.tau(no_cell)     = 0;  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    case 5  % Replication of DIP infected cells
        no_cell = no_cu + no_cv + no_cd;
        no_cell = no_cell +1;  % No of cells increase by 1
        
        % Here we consider that No. of DIPs inside the old cell divides
        % into halves and goes to both cells.
        
        % New cell profile   
        cell.nvir(no_cell) = floor(0.5* cell.nvir(cellwhich));

        % viruses transfering into new cell
        for i =1: cell.nvir(no_cell)
        % Copying vir structures from old cell to new cell    
        cell.vir_pop(no_cell).vir_in(i).type = cell.vir_pop(cellwhich).vir_in(end).type;
        cell.vir_pop(no_cell).vir_in(i).profile  = cell.vir_pop(cellwhich).vir_in(end).profile;    
        % Deleting transfered virus structure 
        cell.vir_pop(cellwhich).vir_in(end) = [];
        end
        
        cell.tau(no_cell) = 0;  % Cell age of new cell
        
        if cell.nvir(no_cell) > 0  % If DIP present on new cell
            cell.type(no_cell)  = 2;  % CD type cell
            no_cd   = no_cd +1;
        else  % new cell doesn't have any DIPs
            no_cu = no_cu +1;
            cell.type(no_cell)    = 1;  % Normal cell
            cell.vir_pop(no_cell) = struct('vir_in',struct('type',{0},'profile',{0}));
        end
        
        % Parent cell profile
        cell.nvir(cellwhich) = ceil(0.5* cell.nvir(cellwhich));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    case 6  % Death of CV
        no_cv   = no_cv -1;  % Number of CV reduce by 1
        
        % All viruses come outside
        for i=1: cell.nvir(cellwhich)
            vir_out.type(no_vir +i)= cell.vir_pop(cellwhich).vir_in(i).type ;
            vir_out.profile(no_vir +i,:)= cell.vir_pop(cellwhich).vir_in(i).profile ;
            
            if cell.vir_pop(cellwhich).vir_in(i).type==2  % DIP increament in case of DIP release
                no_dip = no_dip +1;
            end
        end
        
        no_vir  = no_vir + cell.nvir(cellwhich); % all viruses of the cell come outside
        
        % Deletion of the deleted/dead cell profile
        cell.type(cellwhich)=[];
        cell.nvir(cellwhich)=[];
        cell.vir_pop(cellwhich)=[];
        cell.tau(cellwhich)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 7  % Virus replication inside the cell
        cell.nvir(cellwhich) = cell.nvir(cellwhich) +1; % Number of virus increased by 1
        
        %Initially we copy the whole virus structure, then change it if
        %mutation occurs (see below)
        cell.vir_pop(cellwhich).vir_in(cell.nvir(cellwhich)) = cell.vir_pop(cellwhich).vir_in(viruswhich); 
        
        % Now we have to check for mutation
        isdip = false;
        temp_vir_prof = cell.vir_pop(cellwhich).vir_in(cell.nvir(cellwhich)).profile; %New virus's profile
        
        vir_len = sys_para.hr_kb * temp_vir_prof' ; % length of the chosen virus
        prob_of_mutation = sys_para.del_prob *vir_len /sum(sys_para.hr_kb); % prob of mutation according to virus length
        
        if length(find(temp_vir_prof==1)) >1  % If the virus has more than 1 hr sites, then mutation will occur
            for j=1: sys_para.no_hr_sites
                % Checking whether the genome should be deleted or not
                if ( temp_vir_prof(j)==1 && rand < prob_of_mutation )
                    cell.vir_pop(cellwhich).vir_in(end).profile(j) = 0;
                    isdip =true; % virus becomes DIP
                end
            end
        end
            
        if isdip == true ;  cell.vir_pop(cellwhich).vir_in(end).type = 2; end % virus is DIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 8   % Virus release from cell
        no_vir = no_vir + 1;  % Increase in No. of extracellular virus
        
        % Copying virus details to outside
        vir_out.type(no_vir)= cell.vir_pop(cellwhich).vir_in(viruswhich).type;
        vir_out.profile(no_vir,:)= cell.vir_pop(cellwhich).vir_in(viruswhich).profile ;
        
        if vir_out.type(no_vir)== 2 % DIP
            no_dip = no_dip +1;
        end
        % Changes inside cell
        cell.nvir(cellwhich) = cell.nvir(cellwhich) -1;  % No. of virus inside the cell reduces
        cell.vir_pop(cellwhich).vir_in(viruswhich) =[];   % deletion of virus inside cell
        
        %Now to check whether any more viruses are there inside the cell or not
        if cell.nvir(cellwhich)==0  % No virus inside
            no_cu = no_cu +1; % CU will increase by 1
            no_cv = no_cv -1; % CV will decrease by 1
            cell.type(cellwhich) = 1; % cell becomes CU (1)
            cell.tau(cellwhich) = 0; 
        else  % Still there are virus(normal or, DIP) inside cell
            isdip = 0; j=1;  % Counters to check whether any normal vir is present inside cell
            while (isdip ==0 && j<= cell.nvir(cellwhich) )
                if cell.vir_pop(cellwhich).vir_in(j).type == 1 % Normal virus
                    isdip = 1;
                end
                j = j+1;
            end
            
            if isdip == 0 % NO normal virus is present inside, only DIP inside
                no_cv = no_cv -1; % CV will decrease by 1
                no_cd = no_cd +1; % CD will increase by 1
                cell.type(cellwhich) = 2; % cell becomes CD (2)
                cell.tau(cellwhich) = 0; 
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end   % end of switch events

% Now to increase the infection time of all virus infected cells (CV)
CV_pos_new = find(cell.type ==3);  %Positions of CV
cell.tau(CV_pos_new) = cell.tau(CV_pos_new) + delta_t ;

return