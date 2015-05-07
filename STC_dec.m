function y_dec = STC_dec_V2(Y, H)
N = 2;
[Nr, Nc] = size(Y);
Nt = size(H,2);
% % % % create real source vec
% % % xs = (1:N)'*[1+1j];
% % % xs_real = zeros(2*N,1);
% % % xs_real(1:2:2*N) = real(xs);
% % % xs_real(2:2:2*N) = 1j*imag(xs);
y_dec = zeros(2*N,1);
% C is the codebook
C = [1+1j, -2+1j*2;2+1j*2, 1-1j]; % Alamouti Code
C_vec = reshape(C, Nc*Nt, 1);
C_vec_real = zeros(2*Nt*Nc,1);
C_vec_real(1:2:2*Nt*Nc) = real(C_vec);
C_vec_real(2:2:2*Nt*Nc) = 1j*imag(C_vec);
% create codebook transform
T = zeros(2*Nt*Nc, 2*N);
for k=1:2*Nt*Nc
    if(isreal(C_vec_real(k))==1)
        T(k, 2*abs(C_vec_real(k))-1) = 1 * sign(C_vec_real(k));
    else
        T(k, 2*abs(C_vec_real(k))) = 1 * sign(C_vec_real(k)/1j);
    end
    
end

for k = 1:Nr % for each receive antenna, perform separate processing
    yr = zeros(2*Nc,1);
    yr(1:2:2*Nc) = real(Y(k,:))';
    yr(2:2:2*Nc) = imag(Y(k,:))';
    % equivalent channel
    h = H(k,:);
    H_tilde = zeros(2, 2 * Nt);
    H_tilde(1, 1:2: 2 * Nt) = real(h);
    H_tilde(2, 2:2: 2 * Nt) = real(h);
    H_tilde(1, 2:2: 2 * Nt) = -imag(h);
    H_tilde(2, 1:2: 2 * Nt) = imag(h);
    H_eq = kron(eye(Nc), H_tilde)*T;
%     sum(abs(h.^2)) % to see if the parallelization is happening correctly
%     H_eq'*H_eq
    y_dec = y_dec + H_eq' * yr;
end
y_dec = y_dec(1:2:2*N) + 1j*y_dec(2:2:2*N);
end
