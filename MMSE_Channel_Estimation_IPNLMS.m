function [CH_OTFS_DD_es,Channel_TD_new_es]=MMSE_Channel_Estimation_IPNLMS(Len_chan,K_ff,M,N,x_CE,Train_length_M,Train_length_N,Rx_sig_DFE)

chan_order_L= Len_chan+1;
chan_order_K = 2*K_ff+1 ;

ii=1;
X=zeros((chan_order_K),(chan_order_L));


alpha_CE =0.9;
kezi_CE =0.01; 
deta_CE =0.01;
stepa_CE = 0.6;  %%%注意和LMS的miu不一样，miu要尽量小
W_CE=zeros((chan_order_L)*(chan_order_K),1);
% W= CH_OTFS_DD_new(1:L+1,1:K+1)+rand(size(W))*0.0001;
% stat_M=0;
% stat_N=8;
% end_M=63;
% end_N=25;

for de=0:Train_length_M-1
    for k=0 : Train_length_N-1 % 这里减K的目的是为了保证下面所有训练数据不超过导频
        
        for delay = 0:chan_order_L-1
            for doopler=0:2*K_ff
                
                if doopler<=K_ff
                    doopler_temp=doopler;
                else
                    doopler_temp= N+doopler-2*K_ff-1;
                end
                
                idx = mod(de-delay,M) ;
                idy =  mod(k-doopler_temp,N);
                
                if doopler_temp<=N/2
                    if de>=delay
                        beta= exp(1i*2*pi*(de-delay)*doopler_temp/N/M);
                    else
                        beta= exp(1i*2*pi*(de-delay)*doopler_temp/N/M)*exp(-1i*2*pi*(mod(k-doopler_temp,N))/N);
                    end
                else
                    if de>=delay
                        beta= exp(1i*2*pi*(de-delay)*(doopler_temp-N)/N/M);
                    else
                        beta= exp(1i*2*pi*(de-delay)*(doopler_temp-N)/N/M)*exp(-1i*2*pi*(mod(k-doopler_temp,N))/N);
                    end
                end
                
                X(doopler+1,delay+1) =x_CE(idy+1,idx+1)*beta;
            end
        end
        x_temp_CE= reshape(X,[],1);
        y_temp_CE = W_CE'*x_temp_CE;  % 这里的W要用共轭，否则不收敛
        %         y1(m+1,n+1) = sum(sum((( X.*CH_OTFS_DD_new(1:L+1,N/2-K+2:N/2+K+2)))));
        yy=Rx_sig_DFE.';
        d_CE = yy(k+1,de+1);
        e_CE(ii) = d_CE - y_temp_CE;
        
        
        Lw_CE=(chan_order_L)*(chan_order_K);
        K1_temp_CE =  (1-alpha_CE)/(2*Lw_CE) ;
        Norm_f1_CE =  sum(sum(abs(W_CE ),1));
        
        K_1_CE=K1_temp_CE+(1+alpha_CE)*abs(W_CE)/(2*Norm_f1_CE+kezi_CE);
        UU_CE=reshape(X,[],1);
        
        W_CE = W_CE +stepa_CE*K_1_CE.*(UU_CE)*e_CE(ii)'/(conj(UU_CE)'*(K_1_CE.*conj(UU_CE))+deta_CE);
        W_recorde =W_CE;
        
        ii=ii+1;
    end
end
% plot(abs(e_CE))
% figure()
% plot(10*log(abs(e_CE).^2))
% W_mean=mean(W_recorde(:,end-N_data:end),2);
W_mean_CE =   W_recorde ;
H_CE=reshape(W_mean_CE',chan_order_K,chan_order_L);
CH_OTFS_DD_es=zeros(N,M);
CH_OTFS_DD_es([0:K_ff N-K_ff:N-1]+1,1:chan_order_L)=H_CE;

% 利用估计的信道，来重构原来的信道
Num_sampling=M*N+Len_chan;
doppler_index=[0:N/2-1, -N/2:-1];
Channel_TD_temp =zeros(N,Num_sampling);
Channel_TD_new_es =zeros(Len_chan+1,Num_sampling);
for delay=1:Len_chan+1
    for nv=1:length(doppler_index)
        Channel_TD_temp((nv),:)=CH_OTFS_DD_es(nv,delay)*exp(1j*2*pi *(-Len_chan:-Len_chan-1+Num_sampling)*doppler_index(nv)/M/N  );
    end
    Channel_TD_new_es(delay,:) =sum(Channel_TD_temp,1);
end
% Channel_TD=zeros(L_CP+1,M*N+L_CP);
% [N_path,Num_sampling]=size(Channel_TD);
% for itao=1:taps
%     Channel_TD(delay_taps(itao)+1,:)=chan_coef(itao)*exp(1j*2*pi *(-L_CP:-L_CP-1+Num_sampling)*Doppler_taps(itao)/M/N );
% end
%ee_TD= Channel_TD_new_es-Channel_TD;   % 重构信道误差
%NMSE_TD=sum(sum(abs(ee_TD.^2)))/sum(sum(abs(Channel_TD.^2)));
 