% defining the governing equation for the HNSW wave, u(1), u(3) ...
% represents the displacement for the first, second ... beads. u(2),
% u(4)... represents the velocity of the first, second ... beads.

function uprime = disp_function(t,u)

global mass k delta num_beads n  precomp gamma g k_barrier forward num_barriers

% ft=1/(sqrt(2*pi)*0.001)*exp(-(t-0.004).^2/0.000002);
ft=exp(-(t-0.001).^2/0.0000001);
gt=2*sin(2*pi*10000*t);
ht=ft.*gt;

uprime = zeros(2*num_beads,1);

uprime(1) = u(2);                   %velocity of the first bead
% uprime(2) = 1/mass(1)*(-k(2)*((abs(delta(2)-(u(3)-u(1)))+(delta(2)-(u(3)-u(1))))/2)^n)+precomp/mass(1)+amp*(chirp((t-total_time*tlow),fmin,total_time*thigh,fmax)).*((t<(thigh*total_time))&(t>(tlow*total_time)));     % acceleration of the first bead


uprime(2) = ...
                    1/mass(1)*(ht+precomp-k(2)*ifpos(delta(2)-(u(3)-u(1)))^n-gamma*ifpos_damping(u(2)-u(4),delta(2)-(u(3)-u(1))));
for ii = 2:(num_beads-1)
    uprime(ii*2-1) = u(ii*2);
    uprime(ii*2) = 1/mass(ii)*...
                            (k(ii)*ifpos(delta(ii)-(u(ii*2-1)-u(ii*2-3)))^n...
                           -k(ii+1)*ifpos(delta(ii+1)-(u(ii*2+1)-u(ii*2-1)))^n...
                           +gamma*ifpos_damping(u(ii*2-2)-u(ii*2),delta(ii)-(u(ii*2-1)-u(ii*2-3)))...
                           -gamma*ifpos_damping(u(ii*2)-u(ii*2+2),delta(ii+1)-(u(ii*2+1)-u(ii*2-1)))...
                       )+g-k_barrier*(ifpos(u(ii*2-1))/sqrt(2))^n*sqrt(2)*((u(ii*2-1)<0)*forward+(u(ii*2-1)>0)*(1-forward))*(ii>=32-floor(num_barriers/2))*(ii<=31+ceil(num_barriers/2)); %<0 forward bias
end

uprime(2*num_beads-1) = u(num_beads*2);         % velocity of the last bead
uprime(num_beads*2) =  0;                                    %acceleration of the last bead, wall