function [Wx,TFx,Ifreq,GD,Rep,Chirprate,q,t,f] = TFM(x,fs,s,tf_type,gamma)
%% ˵��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �����Ȩ����(copyright)��
% ���Ľ�(���䣺baowenjie@sjtu.edu.cn)
% Bao Wenjie(email��baowenjie@sjtu.edu.cn)
% �Ϻ���ͨ��ѧ ��еϵͳ���񶯹����ص�ʵ���ң�����Ŷӣ�
% State Key Laboratory of Mechanical System and Vibration, Shanghai Jiao
% Tong University(Li Fucai Team)
% �Ϻ�������������·800�� ��е�붯������ѧԺB¥405��
% Room 405, Building B, School of Mechanical Engineering, 800 Dongchuan Road, Minhang District, Shanghai, 200240, P.R.China
% ����ʹ����ע������
% Please indicate the source of the program
% �ο�����(paper)��Second-Order Synchroextracting Transform with Application to Fault Diagnosis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������Ŀ���
    if (nargin > 5)
        error('����������࣡');
    elseif(nargin == 4)
        gamma = 0;  
    elseif(nargin == 3)
        gamma = 0; tf_type = 'STFT';
    elseif(nargin == 2 || nargin == 1 || nargin == 0 )
        error('ȱ�����������');
    end
%tf_type���
    if ((~strcmp(tf_type, 'SET1')) && (~strcmp(tf_type, 'SET2')))
        error(['tf_type��֧�֣�',tf_type,'���ͣ�',...
            '��������֧�֣�SET1��SET2']);
    end
%�����ֵ
    Wx=0;TFx=0;Ifreq=0;GD=0;q=0;Rep=0;Chirprate=0;t=0;f=0;
%% Ԥ����
%�ж��Ƿ�Ϊ������������ת��
    [xrow,~] = size(x);
    if (xrow~=1)
        x = x';
    end
%% ModifiedʱƵ����
    %%%%%%%%%%%%%���㹫������%%%%%%%%%%%%%%%%
    N = length(x);
    t = (0:N-1)/fs;
    fftx = fftshift(fft(x));
    delta_w = 2*pi/N;
    if mod(N,2)==0
        w = -pi+delta_w*(0:N-1);
        L = N/2+1;
    else
        w = -pi+delta_w/2+delta_w*(0:N-1);
        L = (N+1)/2;
    end   
    Omega = w*fs;%����Ƶ��תģ��Ƶ��   
    f = (0:L-1)*delta_w*fs;%Ƶ�Ʋ�����rad��
    TFx = zeros(N,L);
    dt = 1/fs;
    df = (f(2)-f(1))/2/pi;
%%%%%%%%%%%%%����STFT%%%%%%%%%%%%%%%%
    if (strcmp(tf_type, 'STFT'))
        gf = @(Omega) sqrt(2*s)*pi^(1/4).*exp(-(s*Omega).^2/2);
        Wx = zeros(N,L);
        for ptr = 1:L
            gh = gf(Omega-f(ptr));
            gh = conj(gh);
            xcpsi = ifft(ifftshift(gh .* fftx));
            Wx(:,ptr) = xcpsi;
        end
        TFx = Wx;
        f = f/2/pi;
%%%%%%%%%%%%%����SET1%%%%%%%%%%%%%%%%            
    elseif (strcmp(tf_type, 'SET1'))
        %G(w)
        gf = @(Omega) sqrt(2*s)*pi^(1/4).*exp(-(s*Omega).^2/2);
        Wx = zeros(N,L);
        for ptr = 1:L
            gh = gf(Omega-f(ptr));
            gh = conj(gh);
            xcpsi = ifft(ifftshift(gh .* fftx));
            Wx(:,ptr) = xcpsi;
        end
        %wG(w)
        gf = @(Omega) Omega.*sqrt(2*s)*pi^(1/4).*exp(-(s*Omega).^2/2);
        wWx = zeros(N,L);
        for ptr = 1:L
            gh = gf(Omega-f(ptr));
            gh = conj(gh);
            xcpsi = ifft(ifftshift(gh .* fftx));
            wWx(:,ptr) = xcpsi;
        end
        %����˲ʱƵ��
        Ifreq = zeros(N,L);
        for k = 1:N
            for m = 1:L
                Ifreq(k,m) = f(m)+conj(wWx(k,m)./Wx(k,m));
            end
        end
        Ifreq = Ifreq/2/pi;
        Ifreq(abs(Wx)<gamma)=0;
        %����
        f = f/2/pi;
        for k = 1:N
            for m = 1:L
                fre = min(max(1 + round((real(Ifreq(k,m))-f(1))/df),1),L);
                if (fre == m)
                    TFx(k, fre) = Wx(k,m);
                end
            end
        end
%%%%%%%%%%%%%����SET2%%%%%%%%%%%%%%%%            
    elseif (strcmp(tf_type, 'SET2'))
        %g(t)
        gf = @(Omega) sqrt(2*s)*pi^(1/4).*exp(-(s*Omega).^2/2);
        Wx = zeros(N,L);
        for ptr = 1:L
            gh = gf(Omega-f(ptr));
            gh = conj(gh);
            xcpsi = ifft(ifftshift(gh .* fftx));
            Wx(:,ptr) = xcpsi;
        end
        %dg(t)
        gf = @(Omega) (1j*Omega).*sqrt(2*s)*pi^(1/4).*exp(-(s*Omega).^2/2);
        dWx = zeros(N,L);
        for ptr = 1:L
            gh = gf(Omega-f(ptr));
            gh = conj(gh);
            xcpsi = ifft(ifftshift(gh .* fftx));
            dWx(:,ptr) = xcpsi;
        end
        %ddg(t)
        gf = @(Omega) -(Omega.*Omega).*sqrt(2*s)*pi^(1/4).*exp(-(s*Omega).^2/2);
        ddWx = zeros(N,L);
        for ptr = 1:L
            gh = gf(Omega-f(ptr));
            gh = conj(gh);
            xcpsi = ifft(ifftshift(gh .* fftx));
            ddWx(:,ptr) = xcpsi;
        end
        %tg(t)
        gf = @(Omega) -1j*Omega.*s*s*sqrt(2*s)*pi^(1/4).*exp(-(s*Omega).^2/2);
        tWx = zeros(N,L);
        for ptr = 1:L
            gh = gf(Omega-f(ptr));
            gh = conj(gh);
            xcpsi = ifft(ifftshift(gh .* fftx));
            tWx(:,ptr) = xcpsi;
        end
        %tdg(t)
        gf = @(Omega) (s*s*Omega.*Omega-1).*sqrt(2*s)*pi^(1/4).*exp(-(s*Omega).^2/2);
        tdWx = zeros(N,L);
        for ptr = 1:L
            gh = gf(Omega-f(ptr));
            gh = conj(gh);
            xcpsi = ifft(ifftshift(gh .* fftx));
            tdWx(:,ptr) = xcpsi;
        end
        %����˲ʱƵ��
        Denominator = tdWx.*Wx-tWx.*dWx;
        Numerator = ddWx.*tWx-dWx.*tdWx;
        Numerator2 = dWx.*dWx-ddWx.*Wx;
        p = Numerator./Denominator;
        q = Numerator2./Denominator;
        for ptr = 1:L
            p(:,ptr) = p(:,ptr) + 1i*f(ptr);
        end
        Rep = real(p);
        Ifreq = imag(p)/2/pi;
        Ifreq(abs(Denominator)<gamma)=0;
        q(abs(Denominator)<gamma)=0;
        Chirprate = imag(q);
        %����
        f = f/2/pi;
        for k = 1:N
            for m = 1:L
                fre = min(max(1 + round((real(Ifreq(k,m))-f(1))/df),1),L);
                if (fre == m)
                    TFx(k, fre) = Wx(k,m);
                end
            end
        end        
%%%%%%%%%%%%%����%%%%%%%%%%%%%%%%
    else
        error([stft_type,'��֧�֣�',tf_type]);
    end    
end
