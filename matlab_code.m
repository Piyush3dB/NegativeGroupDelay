



% Sample transfer function H:
n=400;
f=logspace(log10(1),log10(1000),n);  % sample with logarithmic frequency scale
w=2*pi*f;
R=3000; L=100; C=0.098e-6; Rd=30000;

H=1 + 1/Rd*1./( (R+i*w*L).^-1 + i*w*C);

% Creat Bode plot of transfer function H:
figure(1); subplot(2,1,1)
semilogx(f,20*log10(abs(H)))
xlabel('Frequency [Hz]','FontName','Times')
ylabel('Magnitude [dB]','FontName','Times')
title('Magnitude Response','FontName','Times','FontSize',12,'FontWeight','bold')
grid on

subplot(2,1,2)
semilogx(f,180/pi*angle(H))
xlabel('Frequency [Hz]','FontName','Times')
ylabel('Phase [Deg]','FontName','Times')
title('Phase Response','FontName','Times','FontSize',12,'FontWeight','bold')
axis([1 1000 -70 70])
grid on

% Prepare design matrix M used by FDLS:
Fs = 4000;  % sampling rate for digital system

M=[f' abs(H)' 180/pi*angle(H)' ones(n,1)];  % Design matrix for FDLS
fid = fopen('dm.dat','w');
fprintf(fid, '%u %u \n', n, Fs);          % write header (number of data points and sampling rate)
fprintf(fid, '%12.8e %12.8e %12.8e %12.8e \n', M');   % write frequency response data
fclose(fid);


% plot impulse response:
tend = 0.2; 
t=0:1/Fs:tend;    % plot impulse response up to time "tend"
g=R/(2*L); 
wr=sqrt(1/(L*C)-R^2/(4*L^2));
h=exp(-g*t)/(Rd*C).*(cos(wr*t)+R/(2*L*wr)*sin(wr*t));  % impulse response of analog filter
hd=filter(B,A,[1 zeros(1,tend*Fs)]);                   % impulse response of discrete filter

figure(9);
hold off
plot(t,hd*Fs,'r.')
hold on
plot(t,h,'b','LineWidth',2)
axis([0 tend -400 450])
xlabel('Time [s]','FontName','Times')
title('Impulse Responses','FontName','Times','FontSize',13,'FontWeight','bold')
legend('Discrete Filter','Analog Filter');
grid on

% compute frequency response of Hd(w):
z=exp(i*w/Fs); 
zM =[ones(1,n); z.^(-1); z.^(-2)];
Hdw=(B'*zM)./(A'*zM);

##################################################################################################
% plot frequency response:
figure(1); subplot(2,1,1)
hold on
semilogx(f,20*log10(abs(Hdw)),'r')
legend('H(j \omega)','Hd(e^{j \omega / Fs })');

subplot(2,1,2)
hold on
semilogx(f,180/pi*angle(Hdw), 'r')

% compare frequence response with 4ms time advance:
Oneto10hz=1:134; g1hz=20*log10(abs(Hdw(1)));
figure(6); subplot(2,1,1)
semilogx(f(Oneto10hz),20*log10(abs(Hdw(Oneto10hz)))-g1hz,'r',f(Oneto10hz),20*log10(abs(exp(2*pi*i*f(Oneto10hz)/Fs*16))),'g')
xlabel('Frequency [Hz]','FontName','Times')
axis([1 10 -1 1])
subplot(2,1,2)
semilogx(f(Oneto10hz),180/pi*angle(Hdw(Oneto10hz)), 'r',f(Oneto10hz),180/pi*angle(exp(2*pi*i*f(Oneto10hz)/Fs*16)),'g')
xlabel('Frequency [Hz]','FontName','Times')

% plot group delay:
figure(2);
dpH = diff(angle(H))./diff(w);
semilogx(f(2:end),-dpH)
hold on
dpHdw = diff(angle(Hdw))./diff(w);
semilogx(f(2:end),-dpHdw,'r')
title('Group Delay','FontName','Times','FontSize',12,'FontWeight','bold');
xlabel('Frequency [Hz]','FontName','Times')
ylabel('Group Delay [s]','FontName','Times')
legend('H(j \omega)','Hd(e^{j \omga / Fs })');
axis([1 1000 -5.e-3 6e-3])
grid on

% run the modulated exponential pulse through filter:
tmin=-0.15; tmax=.15;  % time boundaries
t=tmin:1/Fs:tmax;
tau=25e-3;              % Gaussian pulse width
p=exp(-(t/tau).^2);     % definition of Gaussian pulse
t0 = find(t>=0);       
p2=[exp(-(t(1:t0(1))/tau).^2) exp(-(t(t0(1)+1:end)/(.1*tau)).^2)];  % A pulse consisting of a wide half for negative time and a narrow half for positive time

f1=100;
x1=p.*cos(2*pi*f1*t);   % input signal x1
f2=1;
x2=p.*cos(2*pi*f2*t);   % input signal x2
x3=p2.*cos(2*pi*f2*t);  % input signal x3

socascade = [B' A';B' A';B' A';B' A'];

f1ind = find(f >= f1);
y1=1/abs(Hdw(f1ind(1)))^4*sosfilt(socascade,x1); % output signal y1=Hd[x1]
f2ind = find(f >= f2);
y2=1/abs(Hdw(f2ind(1)))^4*sosfilt(socascade,x2); % output signal y2=Hd[x2]
y3=1/abs(Hdw(f2ind(1)))^4*sosfilt(socascade,x3); % output signal y3=Hd[x3]

tpmin=-0.08; tpmax=0.08;
figure(3);
hold off
plot(t,x1,'b','LineWidth',2)
hold on
plot(t,y1+2.5,'r','LineWidth',2)
plot(t,p,'b-.',t-0.005,-1.1*p+2.5,'r-.')
axis([tpmin tpmax -1.5 4.])
xlabel('Time [s]','FontName','Times')
title('Input / Output Wave Group at f_1=100 Hz','FontName','Times','FontSize',13,'FontWeight','bold')
legend('Input Signal x_1(t)','Output Signal y_1(t)');
grid on

figure(4);
hold off
plot(t,x2,'b','LineWidth',2)
hold on
plot(t,y2+2,'r','LineWidth',2)
axis([tpmin tpmax -.2 3.4])
xlabel('Time [s]','FontName','Times')
title('Input / Output Wave Group at f_2=1 Hz','FontName','Times','FontSize',12,'FontWeight','bold')
legend('Input Signal x_2(t)','Output Signal y_2(t)');
grid on

figure(5);
hold off
plot(t,x3,'b','LineWidth',2)
hold on
plot(t,y3+2,'r','LineWidth',2)
xlabel('Time [s]','FontName','Times')
title('Response to Truncated Pulse','FontName','Times','FontSize',12,'FontWeight','bold')
legend('Input Signal x_3(t)','Output Signal y_3(t)');
axis([tpmin tpmax -.2 3.2])
grid on


% time advance random signal:
x4=bandlimsig(20000,pi*(1-10/4000));  % generate a random signal bandlimited to 10 Hz @ 4000 Hz Fs.
y4=1/abs(Hdw(1))^4*sosfilt(socascade,x4); % output signal y4=Hd[x4]
t4=1/Fs*(0:length(x4)-1);
figure(7);
plot(t4,x4,'b','LineWidth',2)
hold on
plot(t4,y4,'r','LineWidth',2)
xlabel('Time [s]','FontName','Times')
title('Response to Bandlimited Random Signal','FontName','Times','FontSize',12,'FontWeight','bold')
legend('Input Signal x_4(t)','Output Signal y_4(t)');
axis([4 6 -3 3])
grid on


function x = bandlimsig(N,e)
% Generate 2N+1 samples of a (pi-e)-bandlimited random signal.

f = (pi-e)/pi;               % bandlimit fraction
M = ceil(f*N);               % 2*M+1=number of full band samples
r = randn(2*M+1,1);          % white standard normal signal
x = zeros(2*N+1,1);          % initialize with zero-vector

for k=(-N):1:N
   sincvec=sinc(k*(pi-e)/pi-((-M):1:M)');
   x(k+N+1)=dot(r,sincvec);
end
