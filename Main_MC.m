% Code for KMC simulation of Baculovirus infection and Cell culture passaging.
% During passaging Defective Interfering Particles (DIPs) are produced.
% Version 11-4-19
% Developed By ASHOK DAS/J Chakraborty/S Dutta
clear all
close all
% ==== Declare Global variables and constants =====
% Global Variables:
global no_cu no_cd no_cv no_vir no_dip;

% no_cu = Number of uninfected cells; no_cd = Number of DIP infected cells
% no_cv = Number of virus infected cells; no_vir = Number of extracellular virus
% no_dip == Number of extracellular DIP

% ===== Parameters ====================
% MC Simulation parameters
sim_para.sim_vol    = 1*3.73e-7;  % MC simulation box volume (L)   (no_cell = 5.3*10^5 per mL)
sim_para.i_no_vir   = 1*2000; % Initial no. of viruses.
% Initially all viruses are extracellular and no DIP
sim_para.i_no_cell  = 1*200;  % Initial total no of cells.  Initially all cells are uninfected.

% System parameters
sys_para.del_prob    = 0.02; % Probability of mutation
sys_para.type_event  = 8;    % Types of events (Details in event_select_and_time_step.m)
sys_para.no_hr_sites = 8;    % No of hr sites
sys_para.no_passage  = 1;   % No of passage
sys_para.hr_kb = [7.7, 18.6, 44.1, 23, 4.2, 9.5, 13, 13.8]; % Lenght of genomes in kb

% Constants for various rate processes. See model equations.
% The time unit is taken as hour.  (Kinetic parameters)
kin_para.k_div      =  3.8e-2;  % Rate of division of CU and CD; % Fixed using Estimate_k_div.m

kin_para.k_inf      =  1.1144e-11; % Rate of infection of CU,CD,CV
kin_para.k_dth      =  1e-2;  % Rate of cell death
kin_para.k_dth_beta =  1e-2;  % Rate of cell death

kin_para.k_rep      =  0.156;     % Rate of virus replication
kin_para.k_rel      =  5.4*1e-3;  % Rate of virus release

kin_para.tau_d    = 0;    % After tau_d hrs lysis starts for CV
kin_para.tau_vrel = 0;    % After tau_vrel hrs virus release possible

t_psg =  48; % Hours. Passage time.
% ===================================================

%% ===== Defining the extracellular virus structure =======
% This virus structure accounts only the extracellular viruses.
% When the virus enters a cell, it is considered into the cell structure.
% vir_out --> Type (DIP/Normal) + Profile (hr sites present)
% vir_out.type= 1 => Normal virus; vir_out.type= 2 => DIP
% Example: vir_out.profile = [1 1 1 0 1 1 1 0] means 4th and 8th genomes
% are absent. All other genomes are present.

% Initiating the virus profile
vir_out.type    = ones(sim_para.i_no_vir ,1); % Initially no DIPs
vir_out.profile = ones(sim_para.i_no_vir, sys_para.no_hr_sites);    % Initially all genomes are present.

%% Initialization of no. of different entities
no_cu  = sim_para.i_no_cell;   % Initial no of unifected cells (CU)
no_cv  = 0;                    % Initial no of virus infected cells(CV)
no_cd  = 0;                    % Initial no of DIP infected cells(CD)
no_dip = length(find(vir_out.type==2));  % Initial no of dip outside
no_vir = sim_para.i_no_vir;    % Initial no of total viruses outside

%% Defining and initializing the cell structure.
% The following number code for types of cells will be used:
% Uninfected cell (CU) = 1; Only DIP infected cell (CD) = 2
% Virus infected cell (CV) = 3 (Virus infected cells may have DIP)

cell.type = ones(sim_para.i_no_cell,1);   % Initially all are uninfected cells
cell.nvir = zeros(sim_para.i_no_cell,1);  % Initially no virus
cell.vir_pop = repmat(struct('vir_in',struct('type',{0},'profile',{0})),sim_para.i_no_cell,1);
cell.tau  = zeros(sim_para.i_no_cell,1);  % Post infection time.(by wild virus only)
% =============================================================
exp_psg_counter =1;   % Counter for experimental results matching

for i = 1: sys_para.no_passage
    %% ===== Initialization for ith passage =========================
    tic
    time   = 0;  % process time
    %Initialization of no. of different entities for this passage
    no_cu  = length(find(cell.type==1));   % Initial no of unifected cells (CU)
    no_cv  = length(find(cell.type==3));   % Initial no of virus infected cells(CV)
    no_cd  = length(find(cell.type==2));   % Initial no of DIP infected cells(CD)
    no_dip = length(find(vir_out.type==2));  % Initial no of dip outside
    
    EV_count = zeros(sys_para.type_event,1);   % Counting each events
    tot_no_cell_psg_start = no_cu+no_cv+no_cd; % Tot no of cell present at the start of passage
    
    % Print initial values
    fprintf('Psg: %1.0f Start| Tot_Virus_XCell= %1.0f| DIP_Xcell= %1.0f| Normal_Vir_Xcell= %1.0f \n', i, no_vir, no_dip, (no_vir-no_dip))
    fprintf('Psg: %1.0f Start| Time= %1.3f| CPU time= %1.3f| Tot_cell= %1.0f| CU= %1.0f| CV= %1.0f| CD= %1.0f| MOI= %1.3f \n',...
        i, time, toc, tot_no_cell_psg_start , no_cu, no_cv, no_cd, no_vir/tot_no_cell_psg_start )
    
    % Variables for data saving
    t_prev = 0; del_t = 1;  T_out= 0;
    Cell_out = (no_cu+no_cv+no_cd); Dth_ev = 0;
    Cell_cu = no_cu; Cell_cv = no_cv; Cell_cd = no_cd;
    N_wv = no_vir - no_dip; N_dip = no_dip;
    
    %% ============ Infection Starts ===============================
    while (time <= t_psg && (no_cu+ no_cv+ no_cd)> 0)
        %choosing an event and time increament
        [eventwhich, delta_t, rate_per_cv] = event_select_and_time_step(cell, kin_para, sim_para);
        EV_count(eventwhich) = EV_count(eventwhich)+ 1;
        
        % Time increament
        time = time + delta_t;
        
        %selecting cell and virus for the chosen event
        [cellwhich, viruswhich] = cell_vir_select(cell, eventwhich, rate_per_cv, sys_para);
        
        %Event execution
        [cell, vir_out] = event_execution(eventwhich, cellwhich, viruswhich, cell, vir_out, delta_t, sys_para);
        
        %% Data Saving
        if (time - t_prev)>= del_t
            T_out = [T_out, time];  % Time instance
            Cell_out = [Cell_out, (no_cu+ no_cv+ no_cd)]; % Total cell
            Dth_ev = [Dth_ev EV_count(6)];  % No. of cell death events
            Cell_cu = [Cell_cu no_cu];  % Total no. of CU
            Cell_cv = [Cell_cv no_cv];  % Total no. of CV
            Cell_cd = [Cell_cd no_cd];  % Total no. of CD
            N_wv = [N_wv (no_vir-no_dip)]; % Total no. of extracellular wild-virus
            N_dip = [N_dip no_dip];  % Total no. of extracellular DIP
            t_prev = time;
        end
        
    end
    %% final time data update
    T_out = [T_out, time]; Cell_out = [Cell_out, (no_cu+ no_cv+ no_cd)];
    Dth_ev = [Dth_ev EV_count(6)];
    Cell_cu = [Cell_cu no_cu]; Cell_cv = [Cell_cv no_cv]; Cell_cu = [Cell_cu no_cu];
    N_wv = [N_wv (no_vir-no_dip)]; N_dip = [N_dip no_dip];
    
    % ====== Infection process ends =================
    
    %% Cell number plotting for 1st passage
    if i==1
        figure  % Total no of cell
        plot(T_out, Cell_out/sim_para.i_no_cell,'k-','MarkerSize', 10,'Linewidth',1.5)
        hold on
        plot([24 48 72], [1.933 2.291 2.132], 'k*','MarkerSize', 8,'Linewidth',1.5)  % Exp1
        plot([24 48 72], [1.782 1.724 1.609], 'ko','MarkerSize', 8,'Linewidth',1.5)  % Exp2
        plot([24 48 72], [1.323 1.441 1.592], 'ks','MarkerSize', 8,'Linewidth',1.5)  % Exp3
        % title('Number of living cells')
        xlabel('Time (h)','fontsize',16)
        ylabel('Normalized number of living cells','fontsize',16)
        legend({'MC','Experiment 1','Experiment 2','Experiment 3'},...
            'Location','northwest','fontsize',16)
        legend('boxoff')
        box on
        
        figure % No of dead cells
        plot(T_out, Dth_ev/sim_para.i_no_cell,'k-','MarkerSize', 10,'Linewidth',1.5)
        hold on
        plot([24 48 72], [0.223 0.425 0.51],'k*','MarkerSize', 8,'Linewidth',1.5)  %   Exp1
        plot([24 48 72], [0.707 0.741 0.956],'ko','MarkerSize', 8,'Linewidth',1.5)  % Exp2
        plot([24 48 72], [0.275 0.507 0.637],'ks','MarkerSize', 8,'Linewidth',1.5)  % Exp3
        %  title('Number of dead cells')
        xlabel('Time (h)','fontsize',16)
        ylabel('Fraction of dead cells','fontsize',16)
        legend({'MC','Experiment 1','Experiment 2','Experiment 3'},...
            'Location','northwest','fontsize',16)
        legend('boxoff')
        box on
    end
    
    %% ====== Data Extraction and presentation =======
    % Printing results
    fprintf('Psg: %1.0f End| Tot_Virus_XCell= %1.0f| DIP_Xcell= %1.0f| Wild_Vir_Xcell= %1.0f| WT percent= %1.2f\n', i, no_vir, no_dip, (no_vir-no_dip),(no_vir-no_dip)*100/no_vir)
    fprintf('Psg: %1.0f End| Time= %1.1f| CPU time= %1.1f| Tot_cell= %1.0f| CU= %1.0f| CV= %1.0f| CD= %1.0f| Normalized_cell= %1.4f \n',...
        i, time, toc, no_cu+no_cv+no_cd, no_cu, no_cv, no_cd, (no_cu+no_cv+no_cd)/tot_no_cell_psg_start)
    fprintf('================================================================\n')
    
    % Virus distribution wrt length
    vir_len          = vir_out.profile * sys_para.hr_kb';   % Length of viruses
    avg_vir_len      = mean(vir_len);  % Average length of vectors
    [bin, vir_range] = virus_range_new(vir_len);  %% Not considering viruses of length less than 10
    %   vir_range1 = virus_range(vir_len);    %% bin size =10
    
    % Plotting results
    Exp_psg = [1 5 9 12 16 19 25 28 32];  % Experimental passages
    
    if i==Exp_psg(exp_psg_counter)
        figure%(i)
        title(['Passage ',num2str(i)])
        hold on
        plot(bin,vir_range,'o-')
        
        load(['exp_psg',num2str(Exp_psg(exp_psg_counter)),'.mat'])
        %bar(DIPsize(:,1), 0.01*DIPsize(:,2),'y')
        plot(DIPsize(:,1), 0.01*DIPsize(:,2),'*')
        legend('MC','Exp')
        
        exp_psg_counter = exp_psg_counter + 1; % counter increament
    end
    % Save result
    file = ['MC_PSG-', num2str(i)];
    save([file '.mat'], 'time', 'no_vir', 'no_dip', 'no_cu', 'no_cv', 'no_cd',...
        'cell', 'vir_out', 'avg_vir_len', 'vir_len', 'bin', 'vir_range',...
        'exp_psg_counter','EV_count','T_out','Cell_out','Dth_ev','Cell_cu',...
        'Cell_cv','Cell_cd','N_wv','N_dip')
    
    %% ===== Preparing for the next passage =============
    % Draw 10% of the reamaining system and transfer to 90% fresh system
    [cell, vir_out] = new_sample(cell, vir_out, sim_para);
    
    %% ========== Computational time efficiency checking ==========
    % If simulation takes too long, we half the sim box
    if toc > 600 % 10 minutes
        cell_number  = length(cell.type);  % no of cells
        vir_number   = length(vir_out.type); % no of viruses
        
        del_cell_pos = randperm(cell_number, floor(cell_number/2)); % deletion positions
        del_vir_pos  = randperm(vir_number, floor(vir_number/2)); % deletion positions
        
        % deletion of cells and viruses
        cell.type(del_cell_pos)    = [];
        cell.nvir(del_cell_pos)    = [];
        cell.vir_pop(del_cell_pos) = [];
        cell.tau(del_cell_pos)     = [];
        
        vir_out.type(del_vir_pos)      = [];
        vir_out.profile(del_vir_pos,:) = [];
        
        no_vir = no_vir - length(del_vir_pos);  % update number of virus
        % MC Simulation parameters UPDATE
        sim_para.sim_vol    = 0.5 * sim_para.sim_vol;
        sim_para.i_no_vir   = ceil(0.5 * sim_para.i_no_vir);
        sim_para.i_no_cell  = ceil(0.5 * sim_para.i_no_cell);
    end
    
    % If simulation takes short time, we double the sim box
    if (toc < 20 || length(cell.type) <10 || length(vir_out.type) <100)  % 20 second
        % duplication of cells and viruses
        cell.type    = [cell.type; cell.type];
        cell.nvir    = [cell.nvir; cell.nvir];
        cell.vir_pop = [cell.vir_pop; cell.vir_pop];
        cell.tau     = [cell.tau; cell.tau];
        
        vir_out.type    = [vir_out.type; vir_out.type];
        vir_out.profile = [vir_out.profile; vir_out.profile];
        
        no_vir = 2 *no_vir;  % update number of virus
        % MC Simulation parameters UPDATE
        sim_para.sim_vol    = 2 * sim_para.sim_vol;
        sim_para.i_no_vir   = 2 * sim_para.i_no_vir;
        sim_para.i_no_cell  = 2 * sim_para.i_no_cell;
    end
    % ======== Ready for the next passage ==============
end