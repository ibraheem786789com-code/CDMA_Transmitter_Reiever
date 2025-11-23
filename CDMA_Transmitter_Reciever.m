close all; clc; clear;
%% ---------------------------------------------------------------
% 1. READ TWO AUDIO FILES AND DOWNSAMPLE TO 1200 SAMPLES/SEC
% ---------------------------------------------------------------
[x, gs] = audioread('tum_hi_aana_ringtone.wav');
[y, fs] = audioread('chaand_se_parda.wav');
if size(x, 2) > 1, x = mean(x, 2); end
if size(y, 2) > 1, y = mean(y, 2); end
target_fs = 1200;
if gs ~= target_fs, x = resample(x, target_fs, gs); end
if fs ~= target_fs, y = resample(y, target_fs, fs); end

x = x(1:target_fs*5);   % 5 seconds
y = y(1:target_fs*5);

frame_dur = 0.02;                    % 20 ms per frame
Fs = 1200;
samples_per_frame = round(Fs * frame_dur);
num_frames = ceil(length(x) / samples_per_frame);

user1_total = [];
user2_total = [];
pilot_total = [];    % store detected pilot values (per-symbol) across frames

num_bits = 8;
levels = 2^num_bits;
step_size = 2 / (levels - 1);

trellis = poly2trellis(7, [171 133]);
N = 64;
walsh = hadamard(N);
user_id1 = 3;
user_id2 = 4;

code1 = walsh(user_id1, :);
code2 = walsh(user_id2, :);
traceback = 34;

% ---------------- PROPER PILOT CHANNEL CONFIG ----------------
user_id_pilot = 1;            % Walsh code 0 is row 1 in MATLAB's hadamard matrix
code_pilot = walsh(user_id_pilot, :);  % This is Walsh code 0 (all +1s)

% Power allocation (typical CDMA: pilot gets 15-20% of total power)
pilot_power_ratio = 0.2;      % 20% of total power for pilot
user_power_ratio = (1 - pilot_power_ratio) / 2;  % Equal power for two users

% Scale factors to maintain constant total power
pilot_scale = sqrt(pilot_power_ratio);
user1_scale = sqrt(user_power_ratio);
user2_scale = sqrt(user_power_ratio);

fprintf('Power allocation: Pilot=%.1f%%, User1=%.1f%%, User2=%.1f%%\n', ...
    pilot_power_ratio*100, user_power_ratio*100, user_power_ratio*100);
% ------------------------------------------------------------

%% ---------------------------------------------------------------
% 2. CREATE AUDIO DEVICE WRITER FOR REAL-TIME PLAYBACK
% ---------------------------------------------------------------
deviceWriter = audioDeviceWriter('SampleRate', target_fs);

%% ---------------------------------------------------------------
% 3. PROCESS EACH FRAME
% ---------------------------------------------------------------
for f = 1:num_frames
    fprintf('Processing frame %d/%d...\n', f, num_frames);

    % ----- Extract frame samples -----
    idx_start = (f-1)*samples_per_frame + 1;
    idx_end = min(f*samples_per_frame, length(x));

    x_frame = x(idx_start:idx_end);
    y_frame1 = y(idx_start:idx_end);
    y_frame = y_frame1 / max(abs(y_frame1));
    
    % ---- (1) Quantization ----
    xq = round((x_frame + 1) / step_size);
    yq = round((y_frame + 1) / step_size);

    % ---- (2) Binary conversion ----
    tx_bits1 = reshape(de2bi(xq, 8, 'left-msb').', [], 1);
    tx_bits2 = reshape(de2bi(yq, 8, 'left-msb').', [], 1);

    % ---- (3) Conv. encoding + interleaving ----
    encodedBits1 = convenc(tx_bits1, trellis);
    encodedBits2 = convenc(tx_bits2, trellis);

    rng(30); pattern1 = randperm(length(encodedBits1));
    rng(42); pattern2 = randperm(length(encodedBits2));
    interleaved1 = encodedBits1(pattern1);
    interleaved2 = encodedBits2(pattern2);

    % ---- (4) CDMA spreading (users + pilot) ----
    % number of symbols to produce (columns)
    num_sym = max(length(interleaved1), length(interleaved2));

    % pre-allocate composite spread arrays (1 x N*num_sym)
    spread_bits1 = zeros(1, N * num_sym);
    spread_bits2 = zeros(1, N * num_sym);
    pilot_spread  = zeros(1, N * num_sym);

    for i = 1:num_sym
        % user1 with proper power scaling
        if i <= length(interleaved1)
            bit1 = 2*interleaved1(i) - 1;   % map 0->-1, 1->+1
        else
            bit1 = 0; % padding
        end
        spread_bits1((i-1)*N + (1:N)) = user1_scale * bit1 * code1;

        % user2 with proper power scaling
        if i <= length(interleaved2)
            bit2 = 2*interleaved2(i) - 1;
        else
            bit2 = 0;
        end
        spread_bits2((i-1)*N + (1:N)) = user2_scale * bit2 * code2;

        % PILOT: Constant +1 transmitted every symbol with proper power scaling
        % In real IS-95, pilot uses Walsh 0 and is always ON
        pilot_spread((i-1)*N + (1:N)) = pilot_scale * (1 * code_pilot);
    end

    % PROPER CDMA COMBINING: Sum all channels
    tx_chips = spread_bits1 + spread_bits2 + pilot_spread;
    

%% ------------------- PLOTS: CDMA SPREADING -------------------
if f == 1   % plot only for first frame

    figure('Name','CDMA Spreading & Combining','Position',[200 100 1000 600]);

    subplot(3,1,1);
    plot(spread_bits1(1:500));
    title('User 1 Spread Chips (first 500 samples)');
    xlabel('Chip Index'); ylabel('Amplitude'); grid on;

    subplot(3,1,2);
    plot(spread_bits2(1:500));
    title('User 2 Spread Chips (first 500 samples)');
    xlabel('Chip Index'); ylabel('Amplitude'); grid on;

    subplot(3,1,3);
    plot(tx_chips(1:500));
    title('Composite CDMA Chips (User1 + User2 + Pilot)');
    xlabel('Chip Index'); ylabel('Amplitude'); grid on;

end




    % Analyze composite signal statistics
    composite_power = mean(tx_chips.^2);
    user1_power = mean(spread_bits1.^2);
    user2_power = mean(spread_bits2.^2);
    pilot_power = mean(pilot_spread.^2);
    
    if f == 1  % Display power distribution for first frame
        fprintf('Power distribution: Total=%.3f, User1=%.3f, User2=%.3f, Pilot=%.3f\n', ...
            composite_power, user1_power, user2_power, pilot_power);
    end

    % Normalize composite for modulation stage (maintains relative powers)
    rx = (tx_chips / max(abs(tx_chips))).';
  




%% ------------------- RF TRANSMISSION -------------------
fs_RF = 1e5;       % baseband sampling rate for visualization
fc = 2e3;          % RF carrier frequency (simulation)
sps = 8;           % samples per chip

% Upsample to simulate DAC
rx_upsampled = upsample(rx, sps);

t = (0:length(rx_upsampled)-1)/fs_RF;

% Pulse shaping (rectangular)
pulse = ones(1,sps);
rx_filtered = conv(rx_upsampled, pulse, 'same');



%% ------------------- PLOTS: BASEBAND SIGNAL -------------------
if f == 1
    figure('Name','Baseband Pulse-Shaping','Position',[200 100 1000 500]);

    subplot(2,1,1);
    plot(rx_upsampled(1:300));
    title('Upsampled Chips (first 300 samples)');
    xlabel('Sample Index'); grid on;

    subplot(2,1,2);
    plot(rx_filtered(1:300));
    title('Pulse-Shaped Baseband Signal (first 300 samples)');
    xlabel('Sample Index'); grid on;
end

% Split into I/Q (QPSK mapping: consecutive samples)
L = length(rx_filtered);
if mod(L,2)==1
    rx_filtered = [rx_filtered; 1]; % pad to even
end

I_base = rx_filtered(1:2:end);
Q_base = rx_filtered(2:2:end);


%% ------------------- PLOTS: I(t) AND Q(t) -------------------
if f == 1
    figure('Name','Baseband I/Q','Position',[200 100 1000 500]);

    subplot(2,1,1);
    plot(I_base(1:200),'r');
    title('I(t) Baseband (first 200 samples)');
    xlabel('Sample Index'); ylabel('I'); grid on;

    subplot(2,1,2);
    plot(Q_base(1:200),'b');
    title('Q(t) Baseband (first 200 samples)');
    xlabel('Sample Index'); ylabel('Q'); grid on;
end




% RF upconversion
I_RF = I_base .* cos(2*pi*fc*t(1:length(I_base))).';
Q_RF = Q_base .* sin(2*pi*fc*t(1:length(Q_base))).';
RF_signal = I_RF - Q_RF;  % real RF waveform

% Normalize for plotting
RF_signal = RF_signal / max(abs(RF_signal));

%% ------------------- PLOTS: RF WAVEFORM -------------------
if f == 1
    figure('Name','RF Modulated Signal','Position',[200 100 1000 400]);
    plot(RF_signal(1:1000));
    title('RF Waveform (first 1000 samples)');
    xlabel('Sample Index');
    ylabel('Amplitude');
    grid on;
end

%% ------------------- RF RECEPTION -------------------
% Assume ideal channel (no noise)
rx_RF_received = RF_signal;

% Downconvert to baseband
I_rx = 2 * rx_RF_received .* cos(2*pi*fc*t(1:length(rx_RF_received))).';
Q_rx = -2 * rx_RF_received .* sin(2*pi*fc*t(1:length(rx_RF_received))).';

% Low-pass filter (simple moving average) to remove high-frequency components
lp_filter = ones(1,sps)/sps;
I_bb = conv(I_rx, lp_filter, 'same');
Q_bb = conv(Q_rx, lp_filter, 'same');

%% ------------------- PLOTS: RECEIVED BASEBAND I/Q -------------------
if f == 1
    figure('Name','Received I/Q Baseband','Position',[200 100 1000 500]);

    subplot(2,1,1);
    plot(I_bb(1:200),'r');
    title('Recovered I(t) after Downconversion');
    grid on;

    subplot(2,1,2);
    plot(Q_bb(1:200),'b');
    title('Recovered Q(t) after Downconversion');
    grid on;
end





% Reconstruct baseband signal (interleave I/Q)
rx_baseband = zeros(length(I_bb)+length(Q_bb),1);
rx_baseband(1:2:end) = I_bb;
rx_baseband(2:2:end) = Q_bb;

% Downsample back to one sample per chip
rx_downsampled = rx_baseband(1:sps:end);

%% ------------------- PLOTS: RECOVERED CHIPS -------------------
if f == 1
    figure('Name','Recovered Chips','Position',[200 100 1000 400]);
    plot(rx_downsampled(1:500));
    title('Recovered Chips After RF → Baseband → Downsample (first 500)');
    xlabel('Chip Index'); ylabel('Amplitude'); grid on;
end





%% ------------------- CDMA DESPREADING -------------------
rx_chips = rx_downsampled;

rx1 = rx_chips(:).'; % row vector

    % ---- (5) Despread (users + pilot) ----
    rx_chips = reshape(rx1, N, []);    % N x num_sym
    
    % Despread each channel with proper scaling compensation
    decisions = [code1; code2; code_pilot] * rx_chips;  % 3 x num_sym
    
    % Normalize and apply inverse scaling
    decisions_norm = decisions / N;
    detected_user1_raw = decisions_norm(1, :) / user1_scale;
    detected_user2_raw = decisions_norm(2, :) / user2_scale;
    detected_pilot_raw = decisions_norm(3, :) / pilot_scale;
    
%% ------------------- PLOTS: DESPREADING OUTPUT -------------------
if f == 1
    figure('Name','Despreader Output','Position',[200 100 1000 600]);

    subplot(3,1,1);
    stem(detected_user1_raw(1:50));
    title('User1 Correlator Output (first 50)');
    grid on;

    subplot(3,1,2);
    stem(detected_user2_raw(1:50));
    title('User2 Correlator Output (first 50)');
    grid on;

    subplot(3,1,3);
    plot(detected_pilot_raw(1:200));
    title('Pilot Channel Detection (first 200)');
    grid on;
end
    % Binary decisions
    detected_user1 = double(detected_user1_raw > 0);
    detected_user2 = double(detected_user2_raw > 0);
    detected_pilot = detected_pilot_raw;  % keep analog value for pilot

    % ---- (6) Deinterleaving ----
    % Deinterleave only the length consistent with original encoded bits
    m1 = detected_user1(1:length(interleaved1));
    g1 = detected_user2(1:length(interleaved2));

    rng(30);
    pattern1 = randperm(length(m1));
    deinterleaved1 = zeros(1, length(m1));
    deinterleaved1(pattern1) = m1;

    rng(42);
    pattern2 = randperm(length(g1));
    deinterleaved2 = zeros(1, length(g1));
    deinterleaved2(pattern2) = g1;

    % ---- (7) Viterbi decoding ----
    decodedBits1 = vitdec(deinterleaved1, trellis, traceback, 'trunc', 'hard');
    decodedBits2 = vitdec(deinterleaved2, trellis, traceback, 'trunc', 'hard');

    % ---- (8) Reconstruct analog ----
    bytes1 = reshape(decodedBits1(1:floor(length(decodedBits1)/8)*8), 8, []).';
    n1 = bi2de(bytes1, 'left-msb');
    analog1 = (n1 / 128) - 1;

    bytes2 = reshape(decodedBits2(1:floor(length(decodedBits2)/8)*8), 8, []).';
    n2 = bi2de(bytes2, 'left-msb');
    analog2 = (n2 / 128) - 1;
%% ------------------- PLOTS: RECONSTRUCTED AUDIO -------------------
if f == 1
    figure('Name','Reconstructed Analog (First Frame)','Position',[200 100 1000 500]);

    subplot(2,1,1);
    plot(analog1);
    title('User1 Reconstructed Analog (Frame 1)');
    grid on;

    subplot(2,1,2);
    plot(analog2);
    title('User2 Reconstructed Analog (Frame 1)');
    grid on;
end

    % ---- Store pilot reference for this frame ----
    pilot_total = [pilot_total; detected_pilot(:).']; %#ok<AGROW>

    % ---- (9) Store + Real-Time Playback ----
    user1_total = [user1_total; analog1];
    user2_total = [user2_total; analog2];

    % Play only User 1 in real-time
    deviceWriter(analog1);
end

%% ---------------------------------------------------------------
% 4. CLEANUP & SAVE FILES
% ---------------------------------------------------------------
release(deviceWriter);  % release audio device

audiowrite('user1_output.wav', user1_total, Fs);
audiowrite('user2_output.wav', user2_total, Fs);

% Save pilot reference
if ~isempty(pilot_total)
    pilot_vec = pilot_total(:);
    % The pilot should be approximately +1 (with some noise from other channels)
    fprintf('Pilot statistics: Mean=%.3f, Std=%.3f (should be ~1.0 and small std)\n', ...
        mean(pilot_vec), std(pilot_vec));
    
    % Plot pilot detection performance
    figure('Position', [100, 100, 800, 600]);
    subplot(2,1,1);
    plot(pilot_vec(1:min(1000, length(pilot_vec))));
    title('Pilot Channel Detection (First 1000 symbols)');
    ylabel('Amplitude');
    xlabel('Symbol Index');
    grid on;
    ylim([-2 2]);
    
    subplot(2,1,2);
    histogram(pilot_vec, 50);
    title('Pilot Channel Detection Histogram');
    xlabel('Amplitude');
    ylabel('Count');
    grid on;
    
    % Save pilot as audio for inspection
    reps = ceil(length(user1_total) / length(pilot_vec));
    pilot_audio = repmat(pilot_vec, reps, 1);
    pilot_audio = pilot_audio(1:length(user1_total));
    pilot_audio = pilot_audio / max(abs(pilot_audio));
    audiowrite('pilot_output.wav', pilot_audio, Fs);
end

% Plot composite signal constellation
figure('Position', [100, 100, 800, 400]);
subplot(1,2,1);
plot(real(rx(1:1000)), imag(rx(1:1000)), '.');
title('Composite CDMA Signal Constellation');
xlabel('I-component');
ylabel('Q-component');
grid on;
axis equal;

subplot(1,2,2);
histogram(rx, 50);
title('Composite Signal Amplitude Distribution');
xlabel('Amplitude');
ylabel('Count');
grid on;
disp('✅ Real-time CDMA processing complete with proper Pilot channel.');
numErrors = sum(detected_user1(:) ~= interleaved1(:));
 disp(['Total Bit Errors after RF simulation user1: ', num2str(numErrors)]);
 numErrors = sum(detected_user2(:) ~= interleaved2(:));
 disp(['Total Bit Errors after RF simulation user2: ', num2str(numErrors)]);
