clc
clear
close all
rng('default')
%% 发射数据参数设置
M=128;
N=64;
Train_length_M=32;  % 训练符号长度
Train_length_N=16;
len_train_symbol=Train_length_N*Train_length_M;

all_data=M*N ;
len_symbol = all_data-len_train_symbol;  % 数据块长度
Mod = [2, 4, 8];
Rate=[1/2, 2/3];

Len_chan=40;
N_channel=1; % 通道个数
K_ff = 8;
%训练序列 这个数据多数时候可以用PN代替
M_mod=4;
train_bit = randi(M_mod-1, len_train_symbol, 1);
train_data = qammod(train_bit, M_mod);
train_data = reshape(train_data,Train_length_N,Train_length_M);
x_trian=zeros(N, M);
x_trian(1 : Train_length_N, 1 : Train_length_M) = train_data;

num_channel=10;
SNR=10;

filename='ch1.mat';
load(filename);

for ch_idx=1 : num_channel
%     filename='ch.mat';
%     load(filename);
    for SNR_idx=1:length(SNR)
        for mod_idx = 1:length(Mod)
            M_mod = Mod(mod_idx);
            for nRate =1: length(Rate)
                %% 卷积码
                Rate_code = Rate(nRate);    % 编码速率
                len_ori_bit = len_symbol*mod_idx*Rate_code;     % 原始bit
                len_coded_bit = len_symbol*mod_idx;     % 编码后bit

                
                % 交织
                pos = randperm(len_coded_bit);
                inv_pos = zeros(1, len_coded_bit);
                inv_pos(pos) = 1 : len_coded_bit;
                % 编码
                if Rate_code == 1/3
                    trellis = poly2trellis(3, [7,5,4]);
                elseif Rate_code == 1/2
                    trellis = poly2trellis(7, [171 133]);
                elseif Rate_code == 2/3
                    trellis = poly2trellis([5 4], [23 35 0; 0 5 13]);
                end
                data_info_bit = randi([0,1], len_ori_bit, 1);
                
                codedword = convenc(data_info_bit, trellis); 
                conv_interleave=codedword(pos);
                data_info=bi2de(reshape(conv_interleave, [], mod_idx));
                
                data_info = qammod(data_info, M_mod);

                %导频+数据

                data_info_1 = reshape(data_info(1:Train_length_N*(M-Train_length_M)), Train_length_N, []);
                data_info_2 = reshape(data_info((Train_length_N)*(M-Train_length_M)+1 : end), N-Train_length_N, []);

                data_info_pilot = zeros(N, M);
                data_info_pilot(1 : Train_length_N, 1 : Train_length_M) =  train_data;  %导频
                data_info_pilot(1 : Train_length_N, Train_length_M+1 : M) = data_info_1;
                data_info_pilot(Train_length_N+1 : N, :) = data_info_2;

                data_OTFS = OTFS_modulation(N, M, data_info_pilot);

                max_time_slot = all_data+Len_chan;
                Channel_TD = hmat(1:Len_chan+1,1:max_time_slot);
                CH_OTFS_DD_temp = fft(Channel_TD(:, Len_chan+1 : M : Len_chan+all_data), [], 2)/size(Channel_TD(:, Len_chan+1 : M : Len_chan + all_data), 2);
                CH_OTFS_DD = zeros(M,N);
                CH_OTFS_DD(1 : Len_chan+1, :) = CH_OTFS_DD_temp;

                % add one cp
                data_OTFS_CP = [data_OTFS(N*M-Len_chan+1 : N*M); data_OTFS];
                data_OTFS_CP_filter = 0;
                sigma_2 = 10^(-SNR(SNR_idx)/10);

                for itao = 1:Len_chan+1
                data_OTFS_CP_filter_temp=data_OTFS_CP.*Channel_TD(itao,:).';
                data_OTFS_CP_filter = data_OTFS_CP_filter+circshift([data_OTFS_CP_filter_temp;zeros(Len_chan,1)],(itao)-1);
                end

                rng(ch_idx);
                noise = sqrt(sigma_2/2)*(randn(size(data_OTFS_CP_filter)) + 1i*randn(size(data_OTFS_CP_filter)));
                data_OTFS_CP_faded = data_OTFS_CP_filter + noise;
                % discard cp
                data_OTFS_faded = data_OTFS_CP_faded(Len_chan+1 : Len_chan+(N*M));

                %% OTFS demodulation%%%%
                data_faded = OTFS_demodulation(N, M, data_OTFS_faded);
                data_faded = data_faded.';

                %% demodulation

                % known channel
                x_dect_OTFS_MMSE = OTFS_detection_MMSE_System(data_faded, CH_OTFS_DD, sigma_2);
                x_dect_OTFS_MMSE = x_dect_OTFS_MMSE.';

                data_part1 = reshape(x_dect_OTFS_MMSE(1:Train_length_N, Train_length_M+1 : M), 1, []);
                data_part2 = reshape(x_dect_OTFS_MMSE(Train_length_N+1 : N, :), 1, []);
                data_mod= [data_part1, data_part2];
                
                data_demod=qamdemod(data_mod, M_mod);
                data_bits = reshape(de2bi(data_demod, mod_idx), [], 1);
                data_bits = data_bits(inv_pos); % 解交织
                detec_data = vitdec(data_bits, trellis, 34, 'trunc', 'hard');   % 解卷积码

                BER_known(SNR_idx, (mod_idx-1)*length(Rate)+nRate, ch_idx) = BER_Cacula(detec_data, data_info_bit);

                %         % unknown channel
                %         [CH_OTFS_DD_es, Channel_TD_es] = MMSE_Channel_Estimation_IPNLMS(Len_chan, K_ff, M, N, x_trian, Train_length_M, Train_length_N, data_faded);
                %         x_dect_OTFS_MMSE_es = OTFS_detection_MMSE_System(data_faded, CH_OTFS_DD_es.', sigma_2);
                %         x_dect_OTFS_MMSE_es = x_dect_OTFS_MMSE_es.';
            end
        end
    end
end
