% This function computes the time step length (or, Interval of quissance) &
% chooses the next event. All possible 8 events are:
% Ev_1 = infection of CU; Ev_2 = infection of CV; Ev_3 = infection of CD;
% Ev_4 = Reproduction of CU; Ev_5 = Reproduction of CD;
% Ev_6 = Cell death; Ev_7 = Virus replication; Ev_8 = Virus release from cell;

function [eventwhich, delta_t, rate_per_cv] = event_select_and_time_step(cell, kin_para, sim_para)

global no_vir no_cu no_cv no_cd

%% Rates per CV for events 6,7 & 8
rate_per_cv.death   = zeros(1, no_cv);   % Event 6 - Cell death
rate_per_cv.vrepli  = zeros(1, no_cv);   % Event 7 - Virus replication
rate_per_cv.vrel    = zeros(1, no_cv);   % Event 8 - Virus release
%% Initialization of event rates
event_rate = zeros(1,8);

%% calculating the rate of each event
% Rate of infection of uninfected cells (Rate of Event_1)=
% (k_inf* sim_vol)*(no_cu/ sim_vol)*(no_vir/ sim_vol)= k_inf *(no_cu*no_vir)/sim_vol
% SIMILARLY FOR OTHER EVENTS ALSO
event_rate(1) = kin_para.k_inf *(no_cu*no_vir)/(sim_para.sim_vol);   % Rate of infection of uninfected cells
event_rate(2) = kin_para.k_inf *(no_cv*no_vir)/(sim_para.sim_vol);   % Rate of infection of virus infected cells
event_rate(3) = kin_para.k_inf *(no_cd*no_vir)/(sim_para.sim_vol);   % Rate of infection of DIP infected cells
event_rate(4) = kin_para.k_div * no_cu;               % Rate of multiplication of uninfected cells
event_rate(5) = kin_para.k_div * no_cd;               % Rate of multiplication of DIP infected cells


%% Loop to find the rates of events 6, 7 & 8
CV_pos = find(cell.type ==3);  %Positions of CV

for i = 1:no_cv
    rate_per_cv.vrepli(i) = kin_para.k_rep *cell.nvir(CV_pos(i)); %virus replication inside cell
    
    if (cell.tau(CV_pos(i))) >= kin_para.tau_vrel
        rate_per_cv.vrel(i) = kin_para.k_rel...
            *cell.nvir(CV_pos(i)) *(cell.tau(CV_pos(i)) - kin_para.tau_vrel);  %virus release from cell
        
        if cell.tau(CV_pos(i)) >= kin_para.tau_d
            rate_per_cv.death(i) = kin_para.k_dth *exp(kin_para.k_dth_beta*(cell.tau(CV_pos(i))- kin_para.tau_d)); %cell death
        end
    end
end
event_rate (6) = sum(rate_per_cv.death);
event_rate (7) = sum(rate_per_cv.vrepli);
event_rate (8) = sum(rate_per_cv.vrel);

%%  Choosing an event and finding time step length
cum_sum_event_rate = cumsum(event_rate);   %Cumulative sum of event rates
tot_event_rate     = cum_sum_event_rate(end); %Total event rate
random_no = rand * tot_event_rate;
ind = find(cum_sum_event_rate >= random_no,1);
eventwhich = ind; % Selected event

delta_t = -log(rand) / tot_event_rate ; %Time step

return