clc;
clear all;
close all;


Fe = 3250;
Te = 1/Fe;
Nech = 32;

F1 = 512;
F2 = 640;
FMax = 1536;

time = [0:Te:(Nech-1)*Te];
timeDiscrete = [0:1:Nech-1];
frequency = (timeDiscrete/Nech)*Fe;

%Define the signal here (bandlimited or not ...)
signal = cos(2*pi*F1*(time))+cos(2*pi*F2*(time))+cos(2*pi*FMax*(time));

%Compute DFT
% spectrum=zeros(1,Nech);
% for k = timeDiscrete
%     for l = timeDiscrete
%         spectrum(k+1) = spectrum(k+1) + signal(l+1)*exp(-2*pi*j*l*k/Nech);
%     end
% end
spectrum =  signal*exp(-2*pi*j*timeDiscrete'*timeDiscrete/Nech);

%Compute iDFT
% reconstruction=zeros(1,Nech);
% for k = timeDiscrete
%     for l = timeDiscrete
%         reconstruction(k+1) = reconstruction(k+1) + spectrum(l+1)*exp(2*pi*j*l*k/Nech);
%     end
% end
% reconstruction=reconstruction/Nech;
reconstruction = spectrum*exp(2*pi*j*timeDiscrete'*timeDiscrete/Nech)/Nech;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%    Now interpolation will take place   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define the new resolution
Finterp = 6*Fe;
Tinterp = 1/Finterp;
TimeInterp = [0:Tinterp:(Nech-1)*Te];
NechInterp = length(TimeInterp);
if(mod(NechInterp-Nech,2)==1)
    TimeInterp(NechInterp+1) = TimeInterp(NechInterp)+Te;
    NechInterp = length(TimeInterp);
end
TimeInterpDiscrete = [0:NechInterp-1];
fresampled = [0:Finterp/NechInterp:(NechInterp-1)*Finterp/NechInterp];


%Compute original signal value without any interpolation
signalResampled = cos(2*pi*F1*(TimeInterp))+cos(2*pi*F2*(TimeInterp))+cos(2*pi*FMax*(TimeInterp));
spectrumresampled = zeros(1,NechInterp);

%Compute original signal interpolation through matlab resample function
[P,Q] = rat(Finterp/Fe);
interp_matlab = resample(reconstruction,P,Q);

%Compute original signal interpolation by: first padding the fft
padded_spectrum = zeros(1,NechInterp);
Nzeros = NechInterp-Nech;
padded_spectrum = ifftshift([ zeros(1,floor(Nzeros/2)) fftshift(spectrum) zeros(1,floor(Nzeros/2)+rem(Nzeros,2)) ]);
padded_reconstruction = zeros(1,NechInterp);
% Second: computing the iDFT of the padded spectrum
% for k = TimeInterpDiscrete
%     for l = TimeInterpDiscrete
%         padded_reconstruction(k+1) = padded_reconstruction(k+1) + padded_spectrum(l+1)*exp(2*pi*j*l*k/NechInterp);
%     end
% end
% padded_reconstruction=padded_reconstruction/Nech;
padded_reconstruction = padded_spectrum*exp(2*pi*j*TimeInterpDiscrete'*TimeInterpDiscrete/NechInterp)/(Nech);


%Compute original signal interpolation through shannon interpolation method
interp_shannon=zeros(1,NechInterp);
for k = TimeInterpDiscrete
    for l = timeDiscrete
        interp_shannon(k+1) = interp_shannon(k+1) + reconstruction(l+1)*sinc(Fe*(TimeInterp(k+1)-time(l+1)));
    end
end

%Compute original signal interpolation through dirichlet kernel interpolation method
interp_dirichlet=zeros(1,NechInterp);
for k = TimeInterpDiscrete
    for l = timeDiscrete
        if (TimeInterp(k+1) ~= time(l+1))
            x = 2*pi*Fe*(TimeInterp(k+1)-time(l+1))/NechInterp;
            interp_dirichlet(k+1) = interp_dirichlet(k+1) + reconstruction(l+1)*(sin((NechInterp/2 + 0.5)*x)/sin(0.5*x));
        else
            interp_dirichlet(k+1) = NechInterp * reconstruction(l+1);
            break;
        end
    end
end
interp_dirichlet = interp_dirichlet/NechInterp; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%       Let's print out the result       %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
figure(1);

% Ground truth : deterministic signal is recomputed
plot(TimeInterp,signalResampled,'g');
hold on;
% linear interpolation between subsampled points (matlab tracing tool)
plot(time,(reconstruction),'c');
hold on;
% matlab resample command interpolation
plot(TimeInterp,real(interp_matlab(1:NechInterp)),'r');
hold on;
% Padding the spectrum and reconstructing it 
plot(TimeInterp,real(padded_reconstruction),'m');
hold on;
% Shannon interpolation method
plot(TimeInterp,real(interp_shannon),'b');
hold on;
% Dirichlet kernel interpolation method
plot(TimeInterp,real(interp_dirichlet),'k');

xlabel('Time in s','FontSize',16)
ylabel('Signal value (no unit)','FontSize',16)
title('\it{ Various signal reconstruction from fourier transform }','FontSize',16)
legend('True signal', 'Reconstruction with linear interpolation','Reconstruction with matlab resample function', 'Reconstruction with padded spectrum', 'Reconstruction with Shannon interpolation','Reconstruction with Dirichlet interpolation');

figure(2);

% matlab resample command interpolation
plot(TimeInterp,real(interp_matlab(1:NechInterp))-signalResampled,'r');
hold on;
% Padding the spectrum and reconstructing it 
plot(TimeInterp,real(padded_reconstruction)-signalResampled,'m');
hold on;
% Shannon interpolation method
plot(TimeInterp,real(interp_shannon)-signalResampled,'b');
hold on;
% Dirichlet kernel interpolation method
plot(TimeInterp,real(interp_dirichlet)-signalResampled,'k');

xlabel('Time in s','FontSize',16)
ylabel('Error regarding true signal value (no unit)','FontSize',16)
title('\it{ Interpolation Error }','FontSize',16)
legend('Reconstruction with matlab resample function', 'Reconstruction with padded spectrum', 'Reconstruction with Shannon interpolation', 'Reconstruction with Dirichlet interpolation');
