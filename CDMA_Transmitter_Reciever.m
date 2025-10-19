close all; clc; clear;

%% ---------------------------------------------------------------
% 1. READ TWO AUDIO FILES AND DOWNSAMPLE TO 1200 SAMPLES/SEC
% ---------------------------------------------------------------
[x, gs] = audioread('tum_hi_aana_ringtone.wav');
[y, fs] = audioread('voice2.wav');

if size(x, 2) > 1, x = mean(x, 2); end
if size(y, 2) > 1, y = mean(y, 2); end

target_fs = 1200;
if gs ~= target_fs, x = resample(x, target_fs, gs); end
if fs ~= target_fs, y = resample(y, target_fs, fs); end

x = x(1:target_fs*16);   % 15 seconds
y = y(1:target_fs*15);

frame_dur = 0.02;                    % 20 ms per frame
Fs = 1200;                           
samples_per_frame = round(Fs * frame_dur);
num_frames = ceil(length(x) / samples_per_frame);

user1_total = [];
user2_total = [];

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
    y_frame = y(idx_start:idx_end);

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

    % ---- (4) CDMA spreading ----
    spread_bits1 = [];
    spread_bits2 = [];
    for i = 1:length(interleaved1)
        bit1 = 2*interleaved1(i) - 1;
        spread_bits1 = [spread_bits1 bit1 * code1];
    end
    for i = 1:length(interleaved2)
        bit2 = 2*interleaved2(i) - 1;
        spread_bits2 = [spread_bits2 bit2 * code2];
    end

    rx = spread_bits1 + spread_bits2;

    % ---- (5) Despread ----
    rx_chips = reshape(rx, N, []);
    decisions = [code1; code2] * rx_chips;
    decisions_norm = decisions / N;
    detected_bits = double(decisions_norm > 0);

    m = detected_bits(1, :);
    g = detected_bits(2, :);

    % ---- (6) Deinterleaving ----
    rng(30);
    pattern1 = randperm(length(m));
    deinterleaved1 = zeros(1, length(m));
    deinterleaved1(pattern1) = m;

    rng(42);
    pattern2 = randperm(length(g));
    deinterleaved2 = zeros(1, length(g));
    deinterleaved2(pattern2) = g;

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

disp('✅ Real-time CDMA frame processing and playback complete.');
