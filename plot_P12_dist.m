% This program plots the results for Passage 12.

%% Data Loading
load('exp_psg12.mat') % Experiment results
load('MC_PSG-12.mat') % MC results

bins_present_pos = [1 2 4 5 6 7];% Bins present in experiment results

% create new bins
new_vir_range = zeros(6,1);  
new_vir_range(1) = vir_range(1); new_vir_range(2) = vir_range(2)+vir_range(3);
new_vir_range(3) = vir_range(4); new_vir_range(4) = vir_range(5);
new_vir_range(5) = vir_range(6); new_vir_range(6) = vir_range(7);

%% Line plot
figure
plot(DIPsize(:,1), 0.01*DIPsize(:,2),'k*','MarkerSize', 8,'Linewidth',1.5) % exp
hold on
plot(bin(bins_present_pos),new_vir_range,'ko-','MarkerSize', 10,'Linewidth',1.5)  % MC
title(['Passage-',num2str(12)])
xlabel('Length of viruses','fontsize',16)
ylabel('Number fraction','fontsize',16)
xlim([34 134])
ylim([0 0.4])
legend({'Experiment','MC'},'fontsize',16,'Location','northwest')
legend('boxoff')
box on

%%  Bar plot
bar_mat = zeros(length(bin),2);
bar_mat(bins_present_pos,1) = 0.01*DIPsize(end:-1:1,2);
bar_mat(bins_present_pos,2) = new_vir_range;

figure
bar(bin, bar_mat,2.8)
title('Passage-12','fontsize',14)
xlabel('Length of viruses (kb)','fontsize',16)
ylabel('Number fraction','fontsize',16)
xlim([28 140])
legend({'Experiment','MC'},'fontsize',16,'Location','northwest')
legend('boxoff')
box on



