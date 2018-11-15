% TR=1;
% [hrf,p] = spm_hrf(TR);
hrf=TC(:,1);
Y=fft(hrf);
L=length(hrf);
Fs=1/L;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure
subplot(2,2,1)
plot(hrf)
subplot(2,2,2)
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

hrf_s=[zeros(5,1);hrf];
hrf_s=hrf_s(1:300,:);

Es=randn(size(hrf_s));
c=norm(hrf_s,'fro')/(2*norm(Es,'fro'));
hrf_s=hrf_s+c*Es;
       
Ys=fft(hrf_s);
Ls=length(hrf_s);
Fss=1/Ls;
P2s = abs(Ys/Ls);
P1s = P2s(1:Ls/2+1);
P1s(2:end-1) = 2*P1s(2:end-1);

fs = Fss*(0:(Ls/2))/Ls;
subplot(2,2,3)
plot(hrf_s)
subplot(2,2,4)
plot(fs,P1s)
title('Single-Sided Amplitude Spectrum of X shifted')
xlabel('f (Hz)')
ylabel('|P1(f)|')


figure;plot(f/abs(max(f)),P1);hold on; plot(fs/abs(max(fs)),P1s,'r')

