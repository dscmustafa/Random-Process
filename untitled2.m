clear;
clc;

%--------------------> Read the file entered by the user <-----------------

[file,path] = uigetfile(".mat");

fullName = string(path)+string(file);

Input = load(fullName);
% RandomProcess = Input.X;
% t= Input.t;
Input= struct2cell(Input);

RP = Input{1};
time = Input{2};
% 
% Input2 = load("Sample_Process_2022.mat","X","t");
% X=Input2.sc;
% t=Input2.t;


%==========================================================================

%-------------------> Probability of every sample Func <-------------------
RP_Size = size(RP);

FRB = unique(RP,"rows");
FRB_Size = size(FRB);
C=zeros([FRB_Size(1)],1);
for Ni=1:FRB_Size(1)
    for Nj=1:RP_Size(1)
        if(FRB(Ni,:) == RP(Nj,:))
            C(Ni) = C(Ni)+1;
        end
    end
end

C = C/RP_Size(1);
Confirm = int8(sum(C));
fs = 5000;
[pxx,f] = pwelch(RP,100,90,100,fs);
figure
plot(f,10*log10(pxx))

xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
%==========================================================================

%---------------------------> Ensamble Mean <------------------------------
FRB_P = zeros(FRB_Size);
for i=1:FRB_Size(1)
   FRB_P(i,:) = FRB(i,:)*C(i);
end
Ens_Mean = sum(FRB_P);
figure
plot(time,Ens_Mean);

%==========================================================================

%-----------------> Statistical ACF between ith & jth Samples <------------
Ni=5;
Nj=5;
Si = FRB(:,Ni);
Sj = FRB(:,Nj);
SiSj = Si.*Sj;
SiSjP = SiSj.*C;
ACF = sum(SiSjP);
if (Ni==Nj)
    TAP = ACF;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
All_ACF = zeros(FRB_Size(2));
TolaAvgPowA = zeros(1,width(FRB));
for i=1:FRB_Size(2)
    for j=1:FRB_Size(2)
        Sii = FRB(:,i);
        Sjj = FRB(:,j);
        SiiSjj = Sii.*Sjj;
        SiiSjjP=SiiSjj.*C;
        All_ACF(i,j) = sum(SiiSjjP);
        if(i==j)
            TolaAvgPowA(i) = All_ACF(i,j);
        end
    end
end

figure
mesh(abs(fft(All_ACF)));

figure
plot(TolaAvgPowA);


% =========================================================================

%-----------------> Time Mean of the nth sample func <------------
N = 10;
SampleFunc_N = FRB(N,:);
TimeMean = sum(SampleFunc_N)/(time(FRB_Size(2))-time(1));

% =========================================================================

% % =========================================================================
% 
% %-----------------> Time ACF of the nth sample func <------------
NS = 15;
SampleFunction = FRB(N,:);

Tau = 3;

SFS = width(SampleFunction)-3;
d = zeros(1,SFS);
    for j=1:SFS
        d(j) = SampleFunction(j)*SampleFunction(j+Tau);
    end

    d_time = zeros(1,SFS);
    ti = time(1);
    tf=time(width(time));
    step = (tf - ti)/SFS;
    for i=1:SFS
        d_time(i) = i*step;
    end
g = trapz(d_time,d)/10;

TauMax = width(time) - 1;

SampleFunctionA = FRB(5,:);

TACFGEN = zeros(1,TauMax);
for Tau=1:TauMax
    size = width(SampleFunctionA) - Tau;
        TACFA = zeros(1,size);
    for i=1:size
        TACFA(i) = SampleFunctionA(i)*SampleFunctionA(i+Tau);
    end
    Tf = time(width(time));
    Ti = time(1);
    
    d_time = zeros(1,size);
    step = (Tf - Ti)/size;
    for j = 1:size
        d_time(j) = time(j);
    end
    if(width(TACFA)>1)
         integrated = trapz(d_time,TACFA)/(Tf-Ti);
         TACFGEN(Tau) = integrated;
    else
        TACFGEN(Tau) = TACFA;
    end
end
figure
plot(TACFGEN);









% =========================================================================

fftFRB = zeros(FRB_Size);
avg = zeros(FRB_Size);
for k=1:height(FRB)
    fftFRB(k,:) =abs(fft(FRB(k,:)));
end

fftFRB = fftFRB.^2;
for l=1:height(FRB)
  avg(l,:) = fftFRB(l,:).*C(l);  
end


meanfft = sum(avg);
mf = meanfft./Tf;
 figure
 plot(mf);



% t = 0:0.02:2;
% x=-1:0.02:.98;
% Beta = normrnd(0,1,length(x),1);
% 
% sincos = zeros([length(Beta)],[length(t)]);
% 
% for row = 1:length(Beta)
%     for col = 1:length(t)
%         sincos(row,col) = Beta(row)*sin(2*pi*t(col));
%     end
% end
% % Z= RP.*X;
%   save("sincos.mat","sincos","time");

