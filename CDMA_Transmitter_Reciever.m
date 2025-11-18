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
    m = length(rx);

    % Map to bits (QPSK needs 2 bits per symbol)
    % We need to properly quantize the composite signal
    bit_stream = zeros(m*2,1);
    for k = 1:m
        chip_value = rx(k);
        % Quantize to 4 levels for QPSK
        if chip_value < -0.5
            bits = [0 0];
        elseif chip_value < 0
            bits = [0 1];
        elseif chip_value < 0.5
            bits = [1 0];
        else
            bits = [1 1];
        end
        bit_stream(2*k-1:2*k) = bits;
    end

    %% QPSK Modulation
    hMod = comm.QPSKModulator('BitInput', true, 'PhaseOffset', pi/4);
    tx = step(hMod, bit_stream);

    %% No noise (perfect channel)
    rx3 = tx;

    %% QPSK Demodulation
    hDemod = comm.QPSKDemodulator('BitOutput', true, 'PhaseOffset', pi/4);
    rx_bits = step(hDemod, rx3);

    %% Convert bits back to chip values
    rx_data = zeros(m,1);
    for k = 1:m
        bits = rx_bits(2*k-1:2*k).';
        if isequal(bits, [0 0])
            rx_data(k) = -0.75;  % approximate -1
        elseif isequal(bits, [0 1])
            rx_data(k) = -0.25;  % approximate -0.33
        elseif isequal(bits, [1 0])
            rx_data(k) = 0.25;   % approximate 0.33
        else % [1 1]
            rx_data(k) = 0.75;   % approximate 1
        end
    end
    rx = rx_data;

    % ---- (5) Despread (users + pilot) ----
    rx_chips = reshape(rx, N, []);    % N x num_sym
    
    % Despread each channel with proper scaling compensation
    decisions = [code1; code2; code_pilot] * rx_chips;  % 3 x num_sym
    
    % Normalize and apply inverse scaling
    decisions_norm = decisions / N;
    detected_user1_raw = decisions_norm(1, :) / user1_scale;
    detected_user2_raw = decisions_norm(2, :) / user2_scale;
    detected_pilot_raw = decisions_norm(3, :) / pilot_scale;
    
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
