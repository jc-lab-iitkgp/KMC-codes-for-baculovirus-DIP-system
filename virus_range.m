% This function distributes extracellular viruses into bins (We can
% consider different bin lengths also)

function vir_range = virus_range(vir_len)

global no_vir
vir_range = zeros(14,1); % Initialization

vir_range(1) = length(find(vir_len>= 0 & vir_len <10))/no_vir;
vir_range(2) = length(find(vir_len>= 10 & vir_len <20))/no_vir;
vir_range(3) = length(find(vir_len>= 20 & vir_len <30))/no_vir;
vir_range(4) = length(find(vir_len>= 30 & vir_len <40))/no_vir;
vir_range(5) = length(find(vir_len>= 40 & vir_len <50))/no_vir;
vir_range(6) = length(find(vir_len>= 50 & vir_len <60))/no_vir;
vir_range(7) = length(find(vir_len>= 60 & vir_len <70))/no_vir;
vir_range(8) = length(find(vir_len>= 70 & vir_len <80))/no_vir;
vir_range(9) = length(find(vir_len>= 80 & vir_len <90))/no_vir;
vir_range(10) = length(find(vir_len>= 90 & vir_len <100))/no_vir;
vir_range(11) = length(find(vir_len>= 100 & vir_len <110))/no_vir;
vir_range(12) = length(find(vir_len>= 110 & vir_len <120))/no_vir;
vir_range(13) = length(find(vir_len>= 120 & vir_len <=130))/no_vir;
vir_range(14) = length(find(vir_len > 130))/no_vir;

%plot
%x=[5:10:125 130];
% x=10:10:130;
% figure
% bar(x,vir_range);
% hold on


return