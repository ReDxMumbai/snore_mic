close all;
clear all;
clc;
a = csvread('C:\Users\diptiy\Desktop\MIT MEDIA LABS\Finaldata\yen_fin_8min.csv');      %read your respective data file store it in a variable 

L =length(a);                                                               
samplerate=190;                                                             
fig=1;                                                                          
for i=1:1:L                                                                   
    b(i)=a(i,1);                                                             
end                          
% sample=5;                                                                     
%     for k=0:1:L-sample-1%(floor(L/sample)-1);
%         for j=1:1:sample
%             c(j)=b(k+j);
%         end
%         maxi=max(c);
%         mini=min(c);
%         range(k+1)=maxi-mini;   
%     end
%     
%     
   m=L;
%    
%    
%    L2=L-sample;
%    sample2=10; 
%     for k=0:1:(floor(L2/sample2)-1);
%         s=0;
%         for j=1:1:sample2;
%             d(j)=range(k*sample2+j);
%             s=s+d(j);
%         end
%         envelope(k+1)=s/sample2;   
%     end 
%    
%     
%     L =length(envelope);
% sample=5;    
%     for k=0:1:L-sample-1
%          s=0;
%         for j=1:1:sample
%             e(j)=envelope(k+j);
%             s=s+e(j);
%         end
%         envelope(k+1)=s/sample;    
%     end
%     
   

%     figure(fig); fig=fig+1;
%     subplot(4,1,2);
%     f=(1:1:L)*sample2/samplerate;
%     plot(f,envelope);
%     title('envelope')
%     subplot(4,1,1);
%     g=(1:1:m)/samplerate;
%     plot(g,b);
%     title('signal')
%     
%     
%     for fg=1:1:L
%         if(envelope(fg)>500)
%             square(fg)=1;
%         else
%             square(fg)=0;
%         end
%     end
%     subplot(4,1,3);
%     plot(f,square);
%     title('square');
%     
%     peaks=[zeros(1,L)];
%     q1=0;
%     q2=0;
%     for i=1:1:L-1
%         if (square(i)==0&&square(i+1)==1)
%             q1=i;
%         else if (square(i)==1&&square(i+1)==0)
%                 q2=i;
%             end
%         end
%         if ((q2-q1)>0)
%             Q=[q1,q2];
%             Q1=mean(Q);
%             peaks(floor(Q1))=1;
%         end
%     end
% subplot(4,1,4);    
% plot(f,peaks);
% title('peaks');
%     

%SIGNALS DONE

% s=1
%  for i=20:1:200
%      if(square(floor(i/10))==1)
%          snore(s)=b(i);
%          s=s+1;
%      end
%  end
%  
%  Fsnore=fft(snore-mean(snore));
%  figure(fig);fig=fig+1;
%  fx=(1:1:length(snore));
%  subplot(2,1,1);
%  plot(fx,snore);
%  subplot(2,1,2);
%  plot(fx*(sample rate/length(snore)),abs(fftshift(Fsnore)));
%  title('fft of snore');
 
 
Fb=fft(b);
LFb=length(Fb);
Fsb=fftshift(Fb);
Fsb((LFb/2)-20*LFb/samplerate:1:(LFb/2)+20*LFb/samplerate)=0;
Fb=ifftshift(Fsb);
maxim=max(abs(Fb));



for p=1:length(Fb/2)                                                         %max ka posi nikala
    if abs(Fb(p))==maxim
        maxipos=p;
        break;
    end
end

Freqmax=maxipos*samplerate/length(Fb);
midpass=Freqmax*2/samplerate;


fpass=[midpass-0.05 midpass+0.05];%# passband
fstop=[midpass-0.1 midpass+0.1]; %# frequency where it rolls off to half power
Rpass=4;%# max permissible ripples in stopband (dB)
Astop=80;%# min 40dB attenuation
n=cheb2ord(fpass,fstop,Rpass,Astop);%# calculate minimum filter order to achieve these design requirements
[b1,a1]=cheby2(n,Astop,fstop);
[z,p,k]=cheby2(n,Astop,fstop);
[s,g]=zp2sos(z,p,k);%# create second order sections
% figure(fig);fig=fig+1;
% freqz(b1,a1);
Hd=dfilt.df2sos(s,g);%# create a dfilt object.
filsnore=filter(Hd,b);
%figure(fig);fig=fig+1;
xf=1:1:m;
%plot(xf/samplerate,filsnore);
figure(fig);fig=fig+1;
frefilsnore=fft(filsnore);

subplot(4,1,1);
plot(xf/samplerate,b);
title('signal');

subplot(4,1,2);
plot(xf(2:m/2)*(samplerate/length(filsnore)),(abs(Fb(2:m/2))))

subplot(4,1,3);
plot(xf/samplerate,filsnore);
title('filtered signal');

subplot(4,1,4);
plot(xf(2:m/2)*(samplerate/length(filsnore)),abs(frefilsnore(2:m/2)));
                
%FILTER
                                                                        

sample=5;                                                                     
    for k=0:1:m-sample-1%(floor(L/sample)-1);
        for j=1:1:sample
            c(j)=filsnore(k+j);
        end
        maxi=max(c);
        mini=min(c);
        range(k+1)=maxi-mini;   
    end
    
    
   
   
   
   L2=m-sample;
   sample2=10; 
    for k=0:1:(floor(L2/sample2)-1);
        s=0;
        for j=1:1:sample2;
            d(j)=range(k*sample2+j);
            s=s+d(j);
        end
        envelope(k+1)=s/sample2;   
    end 
   
    
    L =length(envelope);
sample=5;    
    for k=0:1:L-sample-1
         s=0;
        for j=1:1:sample
            e(j)=envelope(k+j);
            s=s+e(j);
        end
        envelope(k+1)=s/sample;    
    end
    
   

    figure(fig); fig=fig+1;
    subplot(3,1,2);
    f=1:1:L;
    plot(f*sample2/samplerate,envelope);
    title('envelope')
    subplot(3,1,1);
    g=1:1:m;
    plot(g/samplerate,filsnore);
    title('signal')
    maxenv=max(envelope);
    
    for fg=1:1:L;
        if(envelope(fg)>0.1*maxenv)
            square(fg)=1;
        else
            square(fg)=0;
        end
    end
    %subplot(4,1,3);
   % plot(f*sample2/samplerate,square);
    %title('square');
    
    peaks=[zeros(1,L)];
    q1=0;
    q2=0;
    xp=1;
    for i=1:1:L-1
         if (square(i)==0&&square(i+1)==1)
            q1=i;
        else if (square(i)==1&&square(i+1)==0)
                q2=i;
            end
        end
        if ((q2-q1)>0)
            Q=[q1,q2];
            Q1=mean(Q);
            peaks(floor(Q1))=1;
            xp=xp+1;
        end
    end
subplot(3,1,3);    
plot(f*sample2/samplerate,peaks);
title('peaks');


b1=1:1:length(peaks);
%plot(b1,peaks);
[y1,x1]=findpeaks(peaks);
start=x1(1);
ends=x1(1);
c1=length(x1);
snorestart=[];
snoreend=[];
for i=1:1:c1-1
    if (((x1(i+1)-x1(i))>2*samplerate))
        snoreend=[snoreend,ends]
        snorestart=[snorestart,start]
        start=x1(i+1);
        
        ends=x1(i+1);
    else if ends==x1(c1-1)
            snoreend=[snoreend,ends]
            snorestart=[snorestart,start]
        else
            ends=x1(i+1);
        end
    end
end




%figure(fig); fig=fig+1;
EnvF=fft(envelope-mean(envelope));
EnvFs=fftshift(abs(EnvF));
%plot(f(2:L/2),EnvF(2:L/2));

lowercutoff=2;
uppercutoff=8;
fpass=[lowercutoff/samplerate uppercutoff/samplerate];%# passband
fstop=[0.95*(lowercutoff/samplerate) (uppercutoff/samplerate)*1.05]; %# frequency where it rolls off to half power
Rpass=1;%# max permissible ripples in stopband (dB)
Astop=90;%# min 40dB attenuation
n=cheb2ord(fpass,fstop,Rpass,Astop);%# calculate minimum filter order to achieve these design requirements
[b1,a1]=cheby2(n,Astop,fstop);
[z,p,k]=cheby2(n,Astop,fstop);
[s,g]=zp2sos(z,p,k);%# create second order sections
% figure(fig);fig=fig+1;
% freqz(b1,a1);
Hd=dfilt.df2sos(s,g);%# create a dfilt object.
filenv=filter(Hd,envelope);
%figure(fig);fig=fig+1;

%plot(xf/samplerate,filsnore);
figure(fig);fig=fig+1;
frefilenv=fft(filenv);

subplot(4,1,1);
plot(f/samplerate,envelope);
title('env');

subplot(4,1,2);
plot(f(2:L/2)*(samplerate/length(filenv)),(abs(EnvF(2:L/2))))

subplot(4,1,3);
plot(f/samplerate,filenv);
title('filtered env');

subplot(4,1,4);
plot(f(2:L/2)*(samplerate/length(filenv)),abs(frefilenv(2:L/2)));


figure(fig); fig=fig+1;
g=1:1:m;
plot(xf/samplerate,filsnore);
title('signal');
hold on
yesno=zeros(1,L);
for i=1:1:length(snoreend)
        yesno(snorestart(i):snoreend(i))=max(abs(filsnore));
end
plot(f*sample2/samplerate,yesno,'color','RED','markerSize',15)
hold off;


%     figure(fig); fig=fig+1;
%     Fb=fft(b-mean(b));
%     Fsb=fftshift(Fb);
%     subplot(3,1,1)
%     plot(g,(abs(Fsb)));
%     title('fft of og signal')
%     
%     N=length(Fb);
%     SFsb=zeros(N,1);
%     unit=samplerate/N;
%     filterA=N/2-floor(1/unit);
%     filterB=N/2-floor(0.1/unit);
%     filterC=N/2+floor(0.1/unit);
%     filterD=N/2+floor(1/unit);
%     for c=filterA:1:filterB
%         SFsb(c)= Fsb(c);
%     end
%     for c=filterC:1:filterD
%         SFsb(c)= Fsb(c);
%     end
%     subplot(3,1,2);
%     SFb=ifftshift(SFsb);
%     plot(g,SFsb);
%     title('filtered');
%     snorefilt=ifft(SFb);
%     subplot(3,1,3);
%     plot(g,snorefilt);
%     title('IFT of filter');
    
    
    
    
    
    
    
    

%     figure(fig); fig=fig+1;
%     F=fft(freq-mean(freq));
%     Fs=fftshift(F);
%     subplot(3,1,1)
%     plot(f,(abs(Fs)));
%     
%     
%      N=length(F);
%     SF=zeros(N,1);
%     unit=(samplerate/sample2)/N;
%     filterA=N/2-floor(1/unit);
%     filterB=N/2+floor(1/unit);
%     for c=filterA:1:filterB
%         SFb(c)= Fs(c);
%     end
%     subplot(3,1,2);
%     SFs= fftshift(SF);
%     plot(f,SFs);
%     title('filtered');
%     snorefilt=ifft(SF);
%     subplot(3,1,3);
%     plot(f,snorefilt);
%     title('IFT of filter');
%     
%     
% Ymax=max(abs(Fs));
% u=1;
% Xmax=[0;0];
% for r=1:1:numel(Fs)
%     if abs(Fs(r))==Ymax
%         Xmax(u)=r
%         u=u+1;
%     end
% end
% frq=((Xmax(2)-Xmax(1))/2)*unit
% 
%     Fil=ifft(F);
%     
%     figure(fig); fig=fig+1;
%     plot(f,(4*abs(Fil)));
%     
%     
    
        
%     Nf =numel(F); xAxisf=linspace(1,10,Nf); widthfilter=4; filterpower=2;
%     filter= exp(-(xAxisf.^2./widthfilter^2).^filterpower);
%     filtertimes=20; F=F.*filter.^filtertimes; subplot(2,1,2)
%     Fs=fftshift(F); plot(f,(abs(Fs)));





