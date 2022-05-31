function [X_detec_MMSE] = OTFS_detection_MMSE_System(Rx_sig, CH_DD, theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x_dect] = OFDM_detection(Rx_sig,CH_FD,OFDM_Parameter)
%
% INPUTS:      Rx_sig: received signal
%              CH_FD: timr invariant channel in frequency
%              OFDM_Parameter: parameter related to OFDM
%
%
% OUTPUT:      x_dect: the detect data symbol.
%
%
% Comments:
%
%
%
% DESCRIPTION: Recover the original data symbol from the received OFDM signal.
%
%
% AUTHOR:           Jianjun Li,
% COPYRIGHT:
% DATE:             06.10.2016
% Last Modified:    06.20.2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Nfft=OFDM_Parameter.Size_of_FFT;
%%Rx_sig----MxN
%%CH_DD-----MxN
%theta------power of noise
%x----- original data 






%% Heff
[Nfft,NumOFDMsubframe]=size(Rx_sig);
M=Nfft;
N=NumOFDMsubframe;

[N_path, N_sampling]=size(CH_DD);
H_DD=zeros(M,N_sampling);
HD_DD_temp=CH_DD;
H_DD(1:N_path,:)=HD_DD_temp;

% Rx_sig_temp=fft(Rx_sig,[],2)*sqrt(NumOFDMsubframe);

y=reshape(Rx_sig,N*M,1);

H_effect1=zeros(N*M,N*M);
H_matrix_all=zeros(M,M,N);  % 一次写出所有的循环矩阵，接下来的矩阵就是在这里面矩阵
for rowIdx=1:N
    temp=H_DD(:,rowIdx);  % 对应的是H_DD中的每个列 也就是对应Doppler
    for j=1:M
        H_matrix_all(:,j,rowIdx)=circshift(temp',[0,j-1])';
    end
end
Location_H=zeros(N,N);  % 上面的子矩阵在总的H中的位置
for j=1:N
    Location_H(:,j )=circshift(1:N, [0,j-1])';
end
for iii=1:N
    index_x=(iii-1)*M+1:iii*M;
    for jjj=1:N
        index_y=(jjj-1)*M+1:jjj*M;
        kkk=Location_H(iii,jjj);
        H_effect1(index_x,index_y)=H_matrix_all(:,:,kkk); % 按照位置填充信道矩阵
    end
end

%% All Phase
All_phase=zeros(M*N,M*N);
delay_orignal=mod(M*1-(M*(1-1):M*1-1),M);

for data_index=1:M*N

    delay=circshift(delay_orignal,data_index-1);  %%
    if mod(data_index,M)==0
        data_delay_index=M;
        data_doopler_index=floor(data_index/M);
        H_index=floor(data_index/M);
    else
        data_delay_index=mod(data_index,M);
        data_doopler_index=floor(data_index/M)+1;
        H_index=floor(data_index/M)+1;
    end
    de=data_delay_index-1;
    phase_division=de-delay;
    phase1_index=find(phase_division>=0); %%de>=delay的相位   0-M
    phase2_index=find(phase_division<0);  %%de<delay的相位    0-M
    for k=1:N
        doopler_block_index=Location_H(   H_index   ,k);
        doopler=doopler_block_index-1;  %%得到块对应的doppler         
        if doopler>N/2
            %de>=delay
            qqqq=phase1_index+M*(k-1);   %All phase中的位置索引
            wwww=phase1_index;                             %delay是循环一维数组，索引只有1：128
            All_phase(data_index,qqqq)  =exp(1i*2*pi*((data_delay_index-1)-delay(wwww))*(doopler-N)/N/M);
            %de<delay
            pppp=phase2_index+M*(k-1);
            oooo=phase2_index;
            All_phase(data_index,pppp)= exp(1i*2*pi*(((data_delay_index-1)-delay(oooo))*(doopler-N)/N/M))*exp(-1i*2*pi*(mod((data_doopler_index-1)-(doopler-N),N))/N);
        else
            %de>=delay
            qqqq=phase1_index+M*(k-1);
            wwww=phase1_index;
            All_phase(data_index,qqqq)  =exp(1i*2*pi*((data_delay_index-1)-delay(wwww))*doopler/N/M);
            %de<delay
            pppp=phase2_index+M*(k-1);
            oooo=phase2_index;
            All_phase(data_index,pppp)= exp(1i*2*pi*(((data_delay_index-1)-delay(oooo))*doopler/N/M))*exp(-1i*2*pi*(mod((data_doopler_index-1)-doopler,N))/N);
        end
    end
end

H_effect=H_effect1.*All_phase;
%Effective data AND Original data
% data_effect=H_effect*xx;  
% data_effect=reshape(data_effect,M,N);
% difference1=data_effect-Rx_sig;
% error_all=sum(sum((abs(difference1.^2))))/sum(sum((abs(Rx_sig.^2))))
% %%在无噪声时，effective data与接受数据的差别----区别Heff是否有错


I=eye(N*M);
W=(H_effect'*H_effect+I*theta)\H_effect';
x_temp=W*y;

clear H_matrix_all All_phase W
clear H_effect H_effect1
X_detec_MMSE=reshape(x_temp,M,N);

