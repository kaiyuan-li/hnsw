% defining the governing equation for the HNSW wave, u(1), u(3) ...
% represents the displacement for the first, second ... beads. u(2),
% u(4)... represents the velocity of the first, second ... beads.

function uprime = hnswDispFun(t,u,param)


% ft=1/(sqrt(2*pi)*0.001)*exp(-(t-0.004).^2/0.000002);
% ft=exp(-(t-0.001).^2/0.0000001);
% gt=2*sin(2*pi*10000*t);
% ht=ft.*gt;

uprime = zeros(2*param.NumOfBeads,1);

uprime(1) = u(2);                   %velocity of the first bead
% uprime(2) = 1/param.MASS(1)*(-k(2)*((abs(param.DELTA(2)-(u(3)-u(1)))+(param.DELTA(2)-(u(3)-u(1))))/2)^n)+param.Precompression/param.MASS(1)+amp*(chirp((t-total_time*tlow),fmin,total_time*thigh,fmax)).*((t<(thigh*total_time))&(t>(tlow*total_time)));     % acceleration of the first bead


uprime(2) = ...
                    1/param.MASS(1)*(param.Precompression-param.K(2)*ifpos(param.DELTA(2)-(u(3)-u(1)))^param.N-param.GAMMA*ifpos_damping(u(2)-u(4),param.DELTA(2)-(u(3)-u(1))));
for idx = 2:(param.NumOfBeads-1)
    uprime(idx*2-1) = u(idx*2);
    uprime(idx*2) = 1/param.MASS(idx)*...
                            (param.K(idx)*ifpos(param.DELTA(idx)-(u(idx*2-1)-u(idx*2-3)))^param.N...
                           -param.K(idx+1)*ifpos(param.DELTA(idx+1)-(u(idx*2+1)-u(idx*2-1)))^param.N...
                           +param.GAMMA*ifpos_damping(u(idx*2-2)-u(idx*2),param.DELTA(idx)-(u(idx*2-1)-u(idx*2-3)))...
                           -param.GAMMA*ifpos_damping(u(idx*2)-u(idx*2+2),param.DELTA(idx+1)-(u(idx*2+1)-u(idx*2-1)))...
                       )+param.G; %<0 forward bias
end

uprime(2*param.NumOfBeads-1) = u(param.NumOfBeads*2);         % velocity of the last bead
uprime(param.NumOfBeads*2) =  0; 