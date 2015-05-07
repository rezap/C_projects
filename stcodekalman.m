clc
clear all
tic

sig2v = 1; % noise variance for all channels
sig2u = 1; % process noise variance
N = 2; % source dimension
Nt = 2; % number of transmitter antennae
Nr = 1; % number of receiver antennae

SNR = 20:2.5:50;
nrep = 1e6;

% % % Eigs = [0.99 0.95]; % eigenvalues of A
% % % D = dctmtx(N);
% % % A = D*diag(Eigs)*D^(-1);
% % % B = sqrt(sig2u)*eye(N);

A = [0.6 -0.8;0.7 0.6];
B = [0.7 0.2;0.2 0.7];


Cu = B*B';

P_th = reshape((eye(N*N) - kron(A,A))^(-1)*reshape(B*B',N*N,1),N,N);%steady-state covariance matrix
p_x_th =  trace(P_th);

d_th = 0.05;
outage = zeros(size(SNR));
outage_th = zeros(size(SNR));
stab_ind = 400;
pow_x = zeros(size(SNR));
%%
parfor l= 1: length(SNR)
    snr = SNR(l);
    gamma = 10^(snr/10);
    x = zeros(N,1) + 1j* zeros(N,1);
    x_hat = zeros(N,1) + 1j* zeros(N,1);
    Mn = ones(N);
    
    
    %     all_d = zeros(1,nrep);
    for k = 1:nrep
        % create source
        u = sqrt(0.5)*(randn(N,1) + 1j*randn(N,1));
        x = A * x + B * u; % x(n) = A x(n-1) + u(n)
        pow_x(l) = pow_x(l) + sum(abs(x.^2));
        [x_STC_enc, Nc] = STC_enc(x); % Nc determines the rate of the code, r = N / Nc
        gamma_sq_const = sqrt(gamma/2/Nt/Nr);
        H = randn(Nr,Nt) + 1j*randn(Nr,Nt); % random fading channel
        %         sum(sum(abs(H.^2)))
        Y =  gamma_sq_const* H * x_STC_enc + sqrt(sig2v/2)*(randn(Nr, Nc) + 1j* randn(Nr, Nc)); % size of Y is Nr x Nc
        y = STC_dec(Y, H); % results in a signal which could have passed through orthogonal channels
        H_norm = sum(sum(abs(H.^2)));
        %%% Perform Kalman filtering
        x_hat = A * x_hat;
        Pn = A * Mn * A' + Cu;
        
        Kn = gamma_sq_const *  Pn * H_norm *(H_norm * sig2v*eye(N) + gamma_sq_const^2 * H_norm^2*Pn)^(-1);
        x_hat = x_hat + Kn * (y - gamma_sq_const * H_norm * x_hat);
        Mn = (eye(N) - Kn * gamma_sq_const * H_norm)*Pn;
        %         all_d(k) = sum(abs((x-x_hat).^2));
        %         dist = sum(abs((x-x_hat).^2));
        dist = 1/N*trace(Mn);
        gamma_dot = H_norm * gamma_sq_const^2;
        % because P(n) is at best equal to C_u, this will result in a lower bound
        dist_th = trace(1/gamma_dot*(eye(2) + 1/gamma_dot * Cu^(-1))^(-1))/N; %
%         dist_th = trace(1/gamma_dot*(eye(2) + 1/gamma_dot * (sig2v*(A*A') + Cu)^(-1))^(-1))/2;
        % This will result in an upper bound
        if ((k > stab_ind) && dist > d_th)
            outage(l) = outage(l) + 1;
        end
        if ((k > stab_ind) && dist_th > d_th)
            outage_th(l) = outage_th(l) + 1;
        end
    end
    outage(l) = outage(l) / (nrep - stab_ind + 1);
    outage_th(l) = outage_th(l) / (nrep - stab_ind + 1);
    
    %     outage(l) = length(find(all_d(stab_ind:end) > d_th)) / (nrep - stab_ind + 1);
end
pow_x = pow_x / (nrep - 1);
time_elapsed = toc;

param.sig2v = sig2v;
param.sig2u = sig2u;
param.N = N;
param.Nt = Nt;
param.Nr = Nr;
param.SNR = SNR;
param.nrep = nrep;
param.Cu = Cu;
param.A = A;
param.B = B;
param.d_th = d_th;
param.outage = outage;
param.outage_th = outage_th;
param.pow_x = pow_x;
param.px_th = p_x_th;
param.time_elapsed = time_elapsed;


semilogy(param.SNR, param.outage,'o-r', 'linewidth',2);
hold on
% semilogy(param.SNR,(param.sig2v*param.Nt*param.Nr*param.px_th/param.d_th)^(param.Nt*param.Nr)...
%     * 10.^(-param.Nt*param.Nr*param.SNR/10),'o-k', 'linewidth',2);
semilogy(param.SNR, param.outage_th,'o-b', 'linewidth',2);
grid on
xlabel('SNR','FontName', 'Times New Roman', 'FontSize', 14)
ylabel('P_{out}','FontName', 'Times New Roman', 'FontSize', 14)
legend('Simulation','Theory')
% % % mu_x = mu_x / nrep;

TIME = clock;
time_str = num2str(TIME(1));
for k = 2:3
    time_str = [time_str '_' num2str(TIME(k))];
end
for k = 4:6
    time_str = [time_str ';' num2str(TIME(k), 2)];
end
% save(['stcodekalman_param_' time_str '.mat'], 'param');


