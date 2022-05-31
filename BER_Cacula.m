%%计算误码率
function [BER]=BER_Cacula(S_code,R_code)
if (length(S_code)~=length(R_code))
    error('The length of S_code and R_code must be same!!!')
else
    BER=length(find(S_code~=R_code))/length(S_code);  %寻找不相等的点
    %     disp(['BER= ',num2str(BER),'%']);
end
% find((S_code-R_code)~=0);