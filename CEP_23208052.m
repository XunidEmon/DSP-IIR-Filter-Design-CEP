% Low-pass IIR Filter Design - CEP EEE 311
% Specs: wp = 0.2*pi, ws = 0.3*pi, Td = 1
% Passband gain: 0.89125 to 1, Stopband gain <= 0.17783

clear; clc; close all;

%% 1. Specifications
Td = 1;
wp = 0.2*pi;          % passband edge (rad/sample)
ws = 0.3*pi;          % stopband edge (rad/sample)

% Convert gain to dB
Gp = -20*log10(0.89125);   % passband ripple in dB (~1 dB)
Gs = -20*log10(0.17783);   % stopband attenuation in dB (~15 dB)

fprintf('Gp = %.4f dB, Gs = %.4f dB\n', Gp, Gs);

%% 2. Bilinear Transform - Pre-warping
Wp = (2/Td)*tan(wp/2);    % pre-warped passband edge
Ws = (2/Td)*tan(ws/2);    % pre-warped stopband edge
fprintf('Wp = %.4f rad/s, Ws = %.4f rad/s\n', Wp, Ws);

%% 3. Butterworth order (Bilinear)
N_bw = ceil(log10((10^(Gs/10)-1)/(10^(Gp/10)-1)) / (2*log10(Ws/Wp)));
fprintf('Butterworth order (Bilinear) N = %d\n', N_bw);

%% 4. Design digital filters (MATLAB built-in)
% --- Butterworth via Bilinear ---
[b_bw, a_bw] = butter(N_bw, wp/pi, 'low');   % 'low' lowpass, cutoff normalized

% --- Chebyshev Type-I via Bilinear ---
N_ch = N_bw;   % same order for fair comparison
Rp = Gp;       % passband ripple in dB
[b_ch, a_ch] = cheby1(N_ch, Rp, wp/pi, 'low');

% --- Chebyshev Type-II via Bilinear ---
Rs = Gs;
[b_ch2, a_ch2] = cheby2(N_ch, Rs, ws/pi, 'low');

%% 5. Impulse Invariance (manual, Butterworth analog)
% Analog Butterworth prototype
[z_a, p_a, k_a] = buttap(N_bw);          % analog prototype
[b_a, a_a] = zp2tf(z_a, p_a, k_a);
[b_a, a_a] = lp2lp(b_a, a_a, Wp);        % scale to Wp
% Convert to digital using impulse invariance
[b_ii, a_ii] = impinvar(b_a, a_a, 1/Td);

%% 6. Frequency response plots
w = linspace(0, pi, 1024);
[H_bw, ~] = freqz(b_bw, a_bw, w);
[H_ch, ~] = freqz(b_ch, a_ch, w);
[H_ch2, ~] = freqz(b_ch2, a_ch2, w);
[H_ii, ~] = freqz(b_ii, a_ii, w);

figure('Name','Magnitude Response Comparison','Color','w');
plot(w/pi, 20*log10(abs(H_bw)), 'b', 'LineWidth', 1.8); hold on;
plot(w/pi, 20*log10(abs(H_ch)), 'r--', 'LineWidth', 1.8);
plot(w/pi, 20*log10(abs(H_ch2)), 'g-.', 'LineWidth', 1.8);
plot(w/pi, 20*log10(abs(H_ii)), 'm:', 'LineWidth', 1.8);
yline(-Gs, 'k--', 'Stopband limit');
xline(0.2, 'k:', 'wp'); xline(0.3, 'k:', 'ws');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');
title('Magnitude Response Comparison of IIR Filters');
legend('Butterworth (Bilinear)','Chebyshev-I','Chebyshev-II','Impulse Invariance', ...
       'Location','southwest');
grid on;

%% 7. Phase response
figure('Name','Phase Response','Color','w');
plot(w/pi, unwrap(angle(H_bw)), 'b', 'LineWidth', 1.8); hold on;
plot(w/pi, unwrap(angle(H_ch)), 'r--', 'LineWidth', 1.8);
plot(w/pi, unwrap(angle(H_ch2)), 'g-.', 'LineWidth', 1.8);
plot(w/pi, unwrap(angle(H_ii)), 'm:', 'LineWidth', 1.8);
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase (rad)');
title('Phase Response Comparison');
legend('Butterworth','Chebyshev-I','Chebyshev-II','Impulse Invariance');
grid on;

%% 8. Pole-Zero Plots
figure('Name','Pole-Zero Plots','Color','w');
subplot(2,2,1);
zplane(b_bw, a_bw); title('Butterworth (Bilinear)');
subplot(2,2,2);
zplane(b_ch, a_ch); title('Chebyshev-I');
subplot(2,2,3);
zplane(b_ch2, a_ch2); title('Chebyshev-II');
subplot(2,2,4);
zplane(b_ii, a_ii); title('Impulse Invariance');

%% 9. Stability check
fprintf('\n--- Stability Check (all |poles| < 1?) ---\n');
fprintf('Butterworth: max pole = %.4f\n', max(abs(roots(a_bw))));
fprintf('Chebyshev-I: max pole = %.4f\n', max(abs(roots(a_ch))));
fprintf('Chebyshev-II: max pole = %.4f\n', max(abs(roots(a_ch2))));
fprintf('Impulse Invariance: max pole = %.4f\n', max(abs(roots(a_ii))));

%% 10. Verify specs
[H_bw_meas, ~] = freqz(b_bw, a_bw, [wp ws]);
fprintf('\nButterworth at wp: %.4f (should be >= 0.89125)\n', abs(H_bw_meas(1)));
fprintf('Butterworth at ws: %.4f (should be <= 0.17783)\n', abs(H_bw_meas(2)));