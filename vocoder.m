% function y=a_loc(x,m1,m2,L,N,TIPO);
m1=0;
m2=14000;
%[x,FS,NBITS]=wavread('lucas.wav');
L=200;
x=cos(2*pi*(1:1000)/25)';
% x=randn(1,1000)';
energia=zeros(floor(length(x)/L),1);
zcr=zeros(floor(length(x)/L)+1,1);
au=zeros(2*L-1,floor(length(x)/L));
ventana=zeros(L,1);
y=[x;zeros((floor(length(x)/L)*L+L)-length(x),1)];

P=10; %orden para coef. Prediccion Lineal
apl=zeros(P,floor(length(x)/L));
%i=ventana tratada
%j=muestra de la ventana i tratada

w=0:4*pi/L:2*pi-4*pi/L;
w=w';


errmatriz=zeros(L,floor(length(x)/L)+1);
auxiliar=0;

% ANALISIS EN VENTANA
%
for i=0:floor(length(x)/L)
    ventana=y(i*L+1:i*L+L);
    
    %Energia, zcr y autocorr
    energia(i+1)=sum(ventana.^2);
    
    for j=2:L
        if sign(2*sign(ventana(j))+1)~=sign(2*sign(ventana(j-1))+1)
            zcr(i+1)=zcr(i+1)+1;
        end
    end

    au(:,i+1)=xcorr(ventana);
    
% CALCULO DFT DE LA VENTANA
%
DFTventana=fft(ventana);
    subplot(2,1,1)
    plot(w/(2*pi)*(FS/2), 20*log(abs(DFTventana(1:length(ventana)/2))))
    subplot(2,1,2)
    plot(ventana), pause



% COEF. PREDICCION LINEAL
%

Rs=xcorr(ventana);
R=Rs(length(ventana):length(ventana)+P-1);
r=Rs(length(ventana)+1:length(ventana)+P);
RT=toeplitz(R);
a=inv(RT)*r;
apl(:,i+1)=a;    
[H,w2]=freqz(1, [1 -a'],512);
%     subplot(2,1,2)
% plot(w2/(2*pi)*(FS/2),20*log(abs(H))), pause

raices=roots([1 -a']);
fases=angle(raices);
    
    
% CODIFICACION DE VOZ
%
err=filter([1 -a'], 1, ventana);
errmatriz(:,i+1)=err;

% subplot(2,1,1)
% plot(errmatriz(:,i+1))
% subplot(2,1,2)
% plot(ventana), pause

end


% for i=0:floor(length(x)/L)   
% err2=filter(1, [1 -apl(:,i+1)']



% energiay2=0;
% 
% for i=1:floor(length(Y)/L)
%     energiay2=[energiay2;energia(i)];
%     energiay2=[energiay2;zeros((L-1),1)];
% end

% plot(y,'r'), hold on
%plot(1:L:length(energia)*L,energia), hold off
% plot(1:L:length(zcr)*L,zcr), hold off
%figure;
%plot(1:L:length(energia)*L,energia+70,'r'), hold on
%plot(1:L:length(zcr)*L,zcr), hold off, grid

%DIBUJAR AUTOCORRELACIONES
% for i=1:floor(length(x)/L)
%     figure;
%     plot(au(:,i)), pause
% end

%DETECTOR ACTIVIDAD VOCAL
actividad=zeros(floor(length(x)/L)+1,1);
for i=1:length(energia)
    if energia(i)>0.1
        actividad(i)=1;
    end
end
sordo=zeros(floor(length(x)/L),1);
for i=1:length(energia)
    if energia(i)>0.08 && energia(i)< 0.15
        sordo(i)=1;
        actividad(i)=0.5;
    end
end



        
% ESTIMACION DEL PITCH
maximos=[0,0];
pitch=0;
pitches=0;
for i=1:length(actividad)
   
    if actividad(i)==1
         maximos=[0,0];
        for j=L+1:length(au(:,i))-1
            if au(j,i)>au(j-1,i) && au(j,i)>au(j+1,i) 
                maximos=[maximos;j,au(j,i)];
            end 
        end
        
        [M,N]=max(maximos(:,2));
        pitch=maximos(N,1)-L;
        pitches=[pitches;pitch];
        
    end
end

pitches2=pitches(2);
for i=3:length(pitches)
    if pitches(i)>50
    pitches2=[pitches2;pitches(i)];
    end
end
estimapitch=FS/(sum(pitches2)/length(pitches2));
pitchmuestras=(sum(pitches2)/length(pitches2))

%DETECCION SONORO/SORDO
%una vez seleccionadas las ventanas donde hay voz(actividad) volvemos a
%comparar con un umbral de energio: niveles bajo->sordos, niveles
%altos->sonoros


%ESPECTROGRAMA
%spectrogram(y,window,0,L,FS)


% CODIFICACION DE VOZ
%
tren=zeros(L,1);
for i=1:floor(pitchmuestras):L
    tren(i)=1;
end
ruido=randn(L,1);

voz_sint=zeros(L,floor(length(x)/L));
for i=0:floor(length(x)/L)
    if actividad(i+1)==1
        voz_sint(:,i+1)=filter(1, [1 -apl(:,i+1)'], tren);
    elseif actividad(i+1)==0.5
        voz_sint(:,i+1)=filter(1, [1 -apl(:,i+1)'], ruido);
    else
        voz_sint(:,i+1)=0.1*ruido;
    end
end

voz_sintetizada=zeros(length(y),1);
for i=0:floor(length(x)/L)
    for j=1:L
        voz_sintetizada(i*L+j)=voz_sint(j,i+1);
    end
end

%wavwrite(y,FS,NBITS,'lucassintetizada.wav');
