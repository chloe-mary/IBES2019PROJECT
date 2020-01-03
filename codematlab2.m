clear all;
close all;

%% 1. Lecture du fichier

[ECGchloegauche, pulsechloegauche] = liremat('IBES/chloe-1left.mat');
[ECGchloedroite, pulsechloedroite] = liremat('IBES/chloe_2right.mat');
[ECGhugogauche, pulsehugogauche] = liremat('IBES/hugo-1left.mat');
[ECGhugodroite, pulsehugodroite] = liremat('IBES/hugo-2right.mat');

%% 2. Selection du signal et affichage
ECG = ECGchloegauche; %qualifier les ondes, leur amplitude etc dans l'article!!
pulse = pulsechloegauche;
fe=1000;
z=(1:(length(ECG)))/fe;
figure(); hold on;
yyaxis left
plot(z, ECG);
ylabel ('Amplitude (in mV)');
yyaxis right
plot(z, pulse);
ylabel ('Amplitude (in mV)');
legend('ECG', 'PPG');
title('ECG and PPG signals without signal processing');
grid;
xlabel('Time (in seconds)');
xlim([0 300]);

%% 3. Filtrage passe-bas 
fcECG=20; %fqce de coupure à 20Hz
fcpulse=20;
ordre = 2;

[b,a]=butter(ordre,fcECG/(fe/2),'low');
freqz(b,a)
ECGfiltered=filter(b,a,ECG);
[b,a]=butter(ordre,fcpulse/(fe/2),'low');
pulsefiltered=filter(b,a,pulse);


figure(); hold on;
yyaxis left
plot(z, ECGfiltered);
ylabel ('Amplitude (in mV)');
yyaxis right
plot(z, pulsefiltered);
ylabel ('Amplitude (in mV)');
legend('ECG', 'PPG');
title('ECG and PPG signal filtered with Butterworth order 2, cut-off frequency of 20Hz');
grid;
xlabel('Time (in seconds)');
xlim([0 300]);

%% 4. Soustraction de la moyenne

largeurdufiltre = 5; %de base 1000
taille=ones(1, largeurdufiltre)/largeurdufiltre;
fDelay = (length(taille)-1)/2;

avglargeurpulsefiltered = filter(taille, 1, pulsefiltered);
pulsefilteredHAUT=pulsefiltered+avglargeurpulsefiltered;
avglargeurECGfiltered = filter(taille, 1, ECGfiltered);
ECGfilteredHAUT=ECGfiltered-avglargeurECGfiltered;
% 
% figure(); hold on;
% yyaxis left
% plot(z, pulsefiltered);
% ylabel ('Amplitude (in mV)');
% yyaxis right
% plot(z, pulsefilteredHAUT);
% plot(z,avglargeurpulsefiltered,'g')
% ylabel ('Amplitude (in mV)');
% legend('PPG', 'PPG centered','PPG average on 5 points');
% title('PPG signal filtered with a Moving Average Filter');
% grid;
% xlabel('Time (in seconds)');
% xlim([0 300]);
% 
% 
% 
% figure(); hold on;
% yyaxis left
% plot(z, ECGfiltered);
% ylabel ('Amplitude (in mV)');
% yyaxis right
% plot(z, ECGfilteredHAUT);
% plot(z,avglargeurECGfiltered,'g')
% ylabel ('Amplitude (in mV)');
% legend('ECG', 'ECG centered','ECG average on 5 points');
% title('ECG signal filtered with a Moving Average Filter');
% grid;
% xlabel('Time (in seconds)');
% xlim([0 300]);
% 
% 
% figure(); hold on;
% yyaxis left
% plot(z, ECGfilteredHAUT);
% ylabel ('Amplitude (in mV)');
% yyaxis right
% plot(z, pulsefilteredHAUT);
% legend('ECG', 'PPG');
% title('Signals ECG and PPG filtered with a Moving Average Filter');
% grid;
% xlabel('Time (in seconds)');
% ylabel ('Amplitude (in mV)');
% xlabel('Time (in seconds)');
% xlim([0 300]);
% 
% 

figure(); hold on;
yyaxis left
plot(ECGfilteredHAUT);
ylabel ('Amplitude (in mV)');
yyaxis right
plot(pulsefilteredHAUT);
legend('ECG', 'PPG');
title('POUR LA SUITE Signals ECG and PPG filtered with a Moving Average Filter');
grid;
xlabel('Time (in seconds)');
ylabel ('Amplitude (in mV)');
xlabel('Samples (in ms)');

%% 4. Extraction des parties exploitables
[x,y] = ginput(2);
ECG=ECGfilteredHAUT(int64(x(1)):int64(x(2)));
pulse=pulsefilteredHAUT(int64(x(1)):int64(x(2)));















% if pulse(1)<ECG(1)
%     pulse(1)=[];
%     ECG(length(ECG))=[];
% end
% while isempty(ECG)
%     [x,y] = ginput(2);
%     ECG=ECGfilteredHAUT(int64(x(1)):int64(x(2)));
%     pulse=pulsefilteredHAUT(int64(x(1)):int64(x(2)));
% end

figure(); hold on;
yyaxis left
plot(ECG);
ylabel ('Amplitude (in mV)');
yyaxis right
plot(pulse);
ylabel ('Amplitude (in mV)');
legend('ECG', 'PPG');
title('ECG and PPG signals filtered and truncated');
grid;
xlabel('Samples (in ms)');


%% 5. Détermination des pics

seuilecg=0.05;
seuilpulse=-1.5;

[pksecg,locsecg] = findpeaks(ECG,'MinPeakDistance',200);
[pkspulse,locspulse] = findpeaks(pulse,'MinPeakDistance',500);

if locspulse(1)<locsecg(1)
    locspulse(1)=[];
    pkspulse(1)=[];
    locsecg(length(locsecg))=[];
    pksecg(length(pksecg))=[];
end

locsecg(pksecg<seuilecg)=[];
pksecg(pksecg<seuilecg)=[];
%locspulse(pkspulse<seuilpulse)=[];
%pkspulse(pkspulse<seuilpulse)=[];

t=1:length(pulse);
x_peaksecg = t(locsecg);
x_peakspulse = t(locspulse);

figure(); hold on;
yyaxis left
plot(t,ECG,x_peaksecg,pksecg,'v')
yyaxis right
plot(t,pulse,x_peakspulse,pkspulse,'v')
legend('ECG', 'pics','pulse', 'pics');
title('ECG and PPG with the peaks succesfully detected');
grid;
xlabel('Samples (in ms)');
ylabel ('Amplitude (in mV)');

%% 6. Calcul du PAT

PATunique = x_peakspulse(3)-x_peaksecg(3) %valeur trop grandes car indices pas adaptés, doit tourner autour de 175

pulsesaute=300; %valeur a régler!!
j=1; %sert d'indice au pulse
i=1; %indice ECG
detection=0; %sert de validation de la détection
PAT=0;
while i<length(x_peaksecg)
    if j>length(x_peakspulse)
        detection=detection+1;
        break
    end
    if x_peakspulse(j)-x_peaksecg(i)<pulsesaute %le pulse a été détecté
        if x_peakspulse(j)-x_peaksecg(i)>0
            PAT=PAT+(x_peakspulse(j)-x_peaksecg(i));
            j=j+1;
            i=i+1;
            detection=detection+1;
        else %pulse détecté avant ECG  %LE 9378 DECONNE!!!
            j=j+1;
            i=i+1;
        end
    else%si pas détecté on incrémente l'ecg pour le comparer avec l'indice du pulse d'aprés
        i=i+1;
    end
end
PATmoyen=PAT/detection