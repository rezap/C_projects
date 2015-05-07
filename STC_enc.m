function [x_STC_enc, Nc] = STC_enc(x)
x_STC_enc = [x(1) -conj(x(2)); x(2) conj(x(1))];
Nc = 2;
end