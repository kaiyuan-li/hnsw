clc
close all
clear all

% initializationclc;close all;clear all;
global mass k delta num_beads n  total_time precomp  gamma g k_barrier forward num_barriers

total_time = 15e-3;         % total time of excitation
num_beads = 200;                  %number of beads, the first one is the striker, the last one is beam or wall
num_barriers = 0;
precomp = 20;                        %precompression force (N)

%steel beads
exx = 195e9*ones(1,num_beads);                % Young's modulus (Pa)
nuxy = 0.3*ones(1,num_beads);                    % Poisson's ratio
dens = 7950*ones(1,num_beads);                    % density (kg/m^3)


initial_velocity=0;

diameter = 19e-3*ones(1,num_beads);        % sphere diameter (m)

exx_barrier = 195e9;
nuxy_barrier = 0.3;

n = 3/2;                            % nonlinear parameter
g = 0;                            % gravity
gamma = 0;        % damping ratio (viscousity wrt speed)

% deducted parameters for the chain of beads
mass = ones(1,num_beads)*4/3*pi.*(diameter/2).^3.*dens;              % mass of each bead
k(1) = inf;                                                          % k(i) is the stiffness between the ith and (i-1)th bead

k(2:num_beads) = 1./(3/2/1.414*((1-nuxy(2:num_beads).^2)./exx(2:num_beads)+(1-nuxy(1:num_beads-1).^2)./exx(1:num_beads-1))).*sqrt(diameter(2:num_beads).*diameter(1:num_beads-1)./(diameter(2:num_beads)+diameter(1:num_beads-1)));

k_barrier = 4*sqrt(diameter(1)/2)/3/((1-nuxy(1)^2)/exx(1)+(1-nuxy_barrier^2)/exx_barrier);
delta = (precomp./k+(1:1:num_beads)*mass(1)*g/k(2)).^(1/n);

% initialize bead related parameters
init_cond = zeros(1,2*num_beads);
init_cond(2) = initial_velocity;


%% forward or reverse
forward = 1; %forward = 1 is forward; =0 is reverse
to_save = 0; %save file to_save = 1; else, not saving


  %% calculate different cases 

    [t,u] = ode15s(@disp_function,(0:5e-7:total_time),init_cond);

%     calculate force between beads
    force=zeros(length(t),num_beads-1);
    for ii = 1:length(t) 
        for jj = 1:num_beads-1
            force(ii,jj) = k(jj+1)*ifpos(delta(jj+1)-(u(ii,jj*2+1)-u(ii,jj*2-1)))^n;
        end
    end

%     calculate force inside beads
    for ii = 1:num_beads-2
        f_bead(:,ii) = (force(:,ii)+force(:,ii+1))/2+eps; 
    end
% forward

fm=max(max(force(:,2:num_beads-2)));
fr=fm/precomp;
v3=0.9314*((4*(exx(1).^2)*precomp)/(diameter(1).^2*dens(1).^3*(1-nuxy(1).^2).^2)).^(1/6)/(fr.^(2/3)-1)*(4/15*(3+2*fr.^(5/3)-5*fr.^(2/3))).^(1/2);
 [m1,t1]=max(f_bead(:,1));
[m2,t2]=max(f_bead(:,60));
 t1=(t1/30001)*total_time;
 t2=(t2/30001)*total_time;
 v4=59*diameter(1)/(t2-t1);
 
 m=0:1:30000;
 N=30001;
 fs=2e6;
 F1 = fft(f_bead(:,1));
 F1 = abs(F1);
 F2 = fft(f_bead(:,num_beads-2));
 F2 = abs(F2);
 f=(m/N)*fs;
 
% draw the colormap for force inside beads
figure(1)
max_force_bead = max(max(f_bead));
max_force_bead1 = max(f_bead(:,1));
max_force_bead2 = max(f_bead(:,num_beads-2));
min_force_bead = min(min(f_bead));
 if (forward==1)
imagesc(1:(num_beads-2), t*1000,(20*log10(abs(f_bead-20)/(max_force_bead-20))))  %fliplr for forward incident
 else
 imagesc(1:(num_beads-2), t*1000,fliplr(abs(f_bead-20)/(max_force_bead-20)))  %fliplr for forward incident
end
axis([0.5 num_beads-1.5 0 total_time*1000])
% colormap(flipud(jet))
colormap(jet);
caxis([-30,0])

colorbar('location','southoutside');
xlabel('Beadsnumber', 'Fontname', 'Times New Roman','FontSize',25)
ylabel('Time(ms)', 'Fontname', 'Times New Roman','FontSize',25)
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 25);
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gcf,'Color',[1 1 1])
 line([31-num_barriers/2 31-num_barriers/2],[0 total_time*1000],'Color','w','LineWidth',2,'linestyle','--')
 line([30+num_barriers/2 30+num_barriers/2],[0 total_time*1000],'Color','w','LineWidth',2,'linestyle','--')

annotation('textbox', [0.15,0.91,0.09,0.1],...
           'String', 'NZ','fontname','times new roman', 'fontsize',25,'linestyle','none');
annotation('textbox', [0.45,0.91,0.11,0.1],...
           'String', 'PBZ','fontname','times new roman', 'fontsize',25,'linestyle','none');
annotation('textbox', [0.77,0.91,0.09,0.1],...
           'String', 'NZ','fontname','times new roman', 'fontsize',25,'linestyle','none');
       annotation('textbox', [0.85,0,0.085,0.1],...
           'String', '(dB)','fontname','times new roman', 'fontsize',25,'linestyle','none');
       

set(gcf,'unit','inch','position',[1 1 6 6]);
set(gca,'YDir','normal','position',[0.135 0.37 0.75 0.57]);      
       
if (to_save==1)
    saveas(gcf,['colormap_reverse_num_barrier_' num2str(num_barriers) '.jpg'])
end

%draw the force-time relationship of a bead
figure(2)
set(gcf, 'units', 'inches', 'position', [6 6 5.5 4])
plot(t*1e3,f_bead(:,1),'k','LineWidth',2) % 31 for 10th last bead
title('The force of the first bead')
xlabel('Time (ms)', 'Fontname', 'Times New Roman','FontSize',22)
ylabel('Force (N)', 'Fontname', 'Times New Roman','FontSize',22)
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 22);
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gcf,'Color',[1 1 1])
% ylim([0 360])
xlim([0,total_time*1000])
if (to_save==1)
 saveas(gcf,['TW_reverse_num_barrier_' num2str(num_barriers) '.emf'])
end
%f_bead60=max(f_bead(:,60));

figure(3)
set(gcf, 'units', 'inches', 'position', [6 6 5.5 4])
plot(t*1e3,f_bead(:,num_beads-2),'k','LineWidth',2) % 31 for 10th last bead
title('The force of the last bead')
xlabel('Time (ms)', 'Fontname', 'Times New Roman','FontSize',22)
ylabel('Force (N)', 'Fontname', 'Times New Roman','FontSize',22)
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 22);
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gcf,'Color',[1 1 1])
% ylim([0 360])
xlim([0,total_time*1000])
if (to_save==1)
 saveas(gcf,['TW_reverse_num_barrier_' num2str(num_barriers) '.emf'])
end

figure(4)
plot(f,F1,'k','LineWidth',2) % 31 for 10th last bead
title('The frequency response of the first bead')
xlabel('Frequency (Hz)', 'Fontname', 'Times New Roman','FontSize',22)
ylabel('Amplitude', 'Fontname', 'Times New Roman','FontSize',22)
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 22);
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gcf,'Color',[1 1 1])
% ylim([0 6e6])
xlim([0,1e4])
if (to_save==1)
 saveas(gcf,['TW_reverse_num_barrier_' num2str(num_barriers) '.emf'])
end

figure(5)
plot(f,F2,'k','LineWidth',2) % 31 for 10th last bead
title('The frequency response of the last bead')
xlabel('Frequency (Hz)', 'Fontname', 'Times New Roman','FontSize',22)
ylabel('Amplitude', 'Fontname', 'Times New Roman','FontSize',22)
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 22);
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gcf,'Color',[1 1 1])
% ylim([0 2.5e6])
xlim([0,6e6])
if (to_save==1)
 saveas(gcf,['TW_reverse_num_barrier_' num2str(num_barriers) '.emf'])
end