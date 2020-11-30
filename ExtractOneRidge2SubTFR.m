function [ExtractTFR,RestTFR,Cs] = ExtractOneRidge2SubTFR(TFx,fs,s,direction,tf_type,delta_or_Rep,q,tradeoff,clear_delta)
%% 说明
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 程序版权归属(copyright)：
% 包文杰(邮箱：baowenjie@sjtu.edu.cn)
% Bao Wenjie(email：baowenjie@sjtu.edu.cn)
% 上海交通大学 机械系统与振动国家重点实验室（李富才团队）
% State Key Laboratory of Mechanical System and Vibration, Shanghai Jiao
% Tong University(Li Fucai Team)
% 上海市闵行区东川路800号 机械与动力工程学院B楼405室
% Room 405, Building B, School of Mechanical Engineering, 800 Dongchuan Road, Minhang District, Shanghai, 200240, P.R.China
% 程序使用请注明出处
% Please indicate the source of the program
% 参考论文(paper)：Second-Order Synchroextracting Transform with Application to Fault Diagnosis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 参数检查
    if (nargin > 9)
        error('输入参数过多！');
    elseif(nargin == 8)
        clear_delta = 0; 
    elseif(nargin == 7)
        clear_delta = 0;         
        tradeoff = 0.009;  
    elseif(nargin == 6 || nargin == 5 || nargin == 4 || nargin == 3 || nargin == 2 || nargin == 1 || nargin == 0)
        error('缺少输入参数！');
    end
%% 预处理
    [N,L] = size(TFx);
    ExtractTFR = zeros(N,L);
    lambda = tradeoff;
    dt = 1/fs;
    df = fs/N;
%% 时间方向压缩的变换重构
    if (strcmp(direction, 'Time') || strcmp(direction, 'T'))
%%%%%%%%%%%%%计算公共部分%%%%%%%%%%%%%%%%
        [Cs,~] = brevridge(abs(TFx), 0:N-1, lambda);
%%%%%%%%%%%%%计算STFT和TSST1和TSST2%%%%%%%%%%%%%%%%
        if (strcmp(tf_type, 'STFT')||strcmp(tf_type, 'TSST1')||strcmp(tf_type, 'TSST2'))
            delta = delta_or_Rep*fs;
            for ptr =1:L
                minvalue = max(Cs(ptr)-round(delta/2),1);
                maxvalue = min(Cs(ptr)+round(delta/2),N);
                ExtractTFR(minvalue:maxvalue,ptr) = TFx(minvalue:maxvalue,ptr);
            end
            RestTFR = TFx - ExtractTFR;
%%%%%%%%%%%%%计算TSET1%%%%%%%%%%%%%%%%            
        elseif (strcmp(tf_type, 'TSET1')||strcmp(tf_type, 'MRM'))
            for ptr =1:L
                ExtractTFR(Cs(ptr),ptr) = TFx(Cs(ptr),ptr);
            end
            RestTFR = TFx - ExtractTFR;
%%%%%%%%%%%%%计算TSET2%%%%%%%%%%%%%%%%            
        elseif (strcmp(tf_type, 'TSET2')||strcmp(tf_type, 'DET'))
            clearTFR = zeros(N,L);
            Rep = delta_or_Rep;
            for ptr =1:L
                ExtractTFR(Cs(ptr),ptr) = TFx(Cs(ptr),ptr);
                minvalue = max(Cs(ptr)-round(clear_delta/2),1);
                maxvalue = min(Cs(ptr)+round(clear_delta/2),N);
                clearTFR(minvalue:maxvalue,ptr) = TFx(minvalue:maxvalue,ptr);
            end
            RestTFR = TFx - clearTFR;
            for ptr =1:L
                M0 = sqrt(-q(Cs(ptr),ptr)+s^2);
                N0 = exp(-0.5*(-Rep(Cs(ptr),ptr)./M0).^2);
                ExtractTFR(Cs(ptr),ptr) = ExtractTFR(Cs(ptr),ptr)*pi^(1/4)/sqrt(s)*M0*N0;
            end
%%%%%%%%%%%%%其它%%%%%%%%%%%%%%%%
        else
            error([direction,'不支持：',tf_type]);
        end
%%%%%%%%%%%%%脊线%%%%%%%%%%%%%%%%
        Cs = (Cs-1)*dt;
%% 频率方向压缩的变换重构       
    else
%%%%%%%%%%%%%计算公共部分%%%%%%%%%%%%%%%%
        [Cs,~] = brevridge(abs(TFx'), 0:L-1, lambda);
%%%%%%%%%%%%%计算STFT和SST1和SST2%%%%%%%%%%%%%%%%
        if (strcmp(tf_type, 'STFT')||strcmp(tf_type, 'SST1')||strcmp(tf_type, 'SST2'))
            delta = delta_or_Rep/fs*N;
            for ptr =1:N
                minvalue = max(Cs(ptr)-round(delta/2),1);
                maxvalue = min(Cs(ptr)+round(delta/2),L);
                ExtractTFR(ptr,minvalue:maxvalue) = TFx(ptr,minvalue:maxvalue);
            end
            RestTFR = TFx - ExtractTFR;
%%%%%%%%%%%%%计算SET1%%%%%%%%%%%%%%%%            
        elseif (strcmp(tf_type, 'SET1')||strcmp(tf_type, 'MRM'))
            for ptr =1:N
                ExtractTFR(ptr,Cs(ptr)) = TFx(ptr,Cs(ptr));
            end
            RestTFR = TFx - ExtractTFR;
%%%%%%%%%%%%%计算SET2%%%%%%%%%%%%%%%%            
        elseif (strcmp(tf_type, 'SET2')||strcmp(tf_type, 'DET'))
            clearTFR = zeros(N,L);
            Rep = delta_or_Rep;
            for ptr =1:N
                ExtractTFR(ptr,Cs(ptr)) = TFx(ptr,Cs(ptr));
                minvalue = max(Cs(ptr)-round(clear_delta/2),1);
                maxvalue = min(Cs(ptr)+round(clear_delta/2),L);
                clearTFR(ptr,minvalue:maxvalue) = TFx(ptr,minvalue:maxvalue);
            end
            RestTFR = TFx - clearTFR;
            for ptr =1:N
                M0 = sqrt(-q(ptr,Cs(ptr))+1/(s^2));
                N0 = exp(-0.5*(-Rep(ptr,Cs(ptr))./M0).^2);
                ExtractTFR(ptr,Cs(ptr)) = ExtractTFR(ptr,Cs(ptr))/pi^(1/4)*sqrt(s)/sqrt(2)*M0*N0;
            end
%%%%%%%%%%%%%其它%%%%%%%%%%%%%%%%
        else
            error([direction,'不支持：',tf_type]);
        end          
    end
%%%%%%%%%%%%%脊线%%%%%%%%%%%%%%%%
    Cs = (Cs-1)*df;
end