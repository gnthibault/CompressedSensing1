clc;
clear all;
close all;



Fe = 3005;
Te = 1/Fe;
Nech = 100;

F1 = 500;
F2 = 1000;
FMax = 1500;

time = [0:Te:(Nech-1)*Te];
timeDiscrete = [1:1:Nech];
frequency = (timeDiscrete/Nech)*Fe;

signal = cos(2*pi*F1*(time))+cos(2*pi*F2*(time))+cos(2*pi*FMax*(time));

%Compute the FFT
spectrum=zeros(1,Nech);
for k = timeDiscrete
    for l = timeDiscrete
        spectrum(k) = spectrum(k) + signal(l)*exp(-2*pi*j*l*k/Nech);
    end
end

%Compute de inverse FFT
reconstruction=zeros(1,Nech);
for k = timeDiscrete
    for l = timeDiscrete
        reconstruction(k) = reconstruction(k) + spectrum(l)*exp(2*pi*j*l*k/Nech);
    end
end
reconstruction=reconstruction/Nech;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%    Now interpolation will take place   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Finterp = 6*Fe;
Tinterp = 1/Finterp;
TimeInterp = [0:Tinterp:(Nech-1)*Te];
[m,n] = size(TimeInterp);
NechInterp = n;
TimeInterpDiscrete = [1:1:NechInterp];

%Compute original signal value without any interpolation
signalResampled = cos(2*pi*F1*(TimeInterp))+cos(2*pi*F2*(TimeInterp))+cos(2*pi*FMax*(TimeInterp));

%Compute original signal interpolation through patlab resample function
[P,Q] = rat(Finterp/Fe);
interp_matlab = resample(reconstruction,P,Q);

%Compute original signal interpolation through shannon interpolation method
interp_shannon=zeros(1,NechInterp);
for k = TimeInterpDiscrete
    for l = timeDiscrete
        interp_shannon(k) = interp_shannon(k) + reconstruction(l)*sinc(Fe*(TimeInterp(k)-time(l)));
    end
end

%Compute original signal interpolation through dirichlet kernel interpolation method
interp_dirichlet=zeros(1,NechInterp);
for k = TimeInterpDiscrete
    for l = timeDiscrete
        if (TimeInterp(k) ~= time(l))
            x = 2*pi*Fe*(TimeInterp(k)-time(l))/NechInterp;
            interp_dirichlet(k) = interp_dirichlet(k) + reconstruction(l)*(sin((NechInterp/2 + 0.5)*x)/sin(0.5*x));
        else
            interp_dirichlet(k) = NechInterp * reconstruction(l);
            break;
        end
    end
end
interp_dirichlet = interp_dirichlet/NechInterp; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%       Let's print out the result       %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(1);
% plot(time,signal);
% hold on;
% plot(time,real(reconstruction),'r');
% 
% figure(2);
% plot(frequency(1:Nech/2),abs(spectrum(1:Nech/2)));


figure(3);

% Ground truth : deterministic signal is recomputed
plot(TimeInterp,signalResampled,'g');
hold on;
% linear interpolation between subsampled points (matlab tracing tool)
plot(time,real(reconstruction),'c');
hold on;
% matlab resample command interpolation
plot(TimeInterp,real(interp_matlab(1:NechInterp)),'r');
hold on;
% Shannon interpolation method
plot(TimeInterp,real(interp_shannon),'b');
hold on;
% Dirichlet kernel interpolation method
plot(TimeInterp,real(interp_dirichlet),'k');

xlabel('Time in s','FontSize',16)
ylabel('Signal value (no unit)','FontSize',16)
title('\it{ Various signal reconstruction from fourier transform }','FontSize',16)
legend('True signal', 'Reconstruction with linear interpolation','Reconstruction with matlab resample function', 'Reconstruction with Shannon interpolation','Reconstruction with Dirichlet interpolation');

figure(4);

% matlab resample command interpolation
plot(TimeInterp,real(interp_matlab(1:NechInterp))-signalResampled,'r');
hold on;
% Shannon interpolation method
plot(TimeInterp,real(interp_shannon)-signalResampled,'b');
hold on;
% Dirichlet kernel interpolation method
plot(TimeInterp,real(interp_dirichlet)-signalResampled,'k');

xlabel('Time in s','FontSize',16)
ylabel('Error regarding true signal value (no unit)','FontSize',16)
title('\it{ Interpolation Error }','FontSize',16)
legend('Reconstruction with matlab resample function', 'Reconstruction with Shannon interpolation', 'Reconstruction with Dirichlet interpolation');
