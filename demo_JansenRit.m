clear all
clf

clc;
%% Runge Kutta's paramters (4th-order)
ra = [0, 0.5, 0.5, 1];
rb = [1, 2, 2, 1];
%%
% beta oscillation
% A = 6.5;
% B = 30;
% C = 135; % # of interneurons
% v0 = 6;
% e0 = 2.5;
% R = 0.56;
% 
% a = 80;
% b = 68;

% % alpha oscillation (Jansen & Rit, 1995)
A = 3.25;
B = 22;
C = 135; % # of interneurons
v0 = 6;
e0 = 2.5;
R = 0.56;

a = 100;
b = 50;

% theta oscillation
% A = 3.25;
% B = 22;
% C = 135; % # of interneurons
% v0 = 6;
% e0 = 2.5;
% R = 0.56;
%
% a = 20;
% b = 17;

C1 = C;
C2 = 0.8*C1;
C3 = C1/4;
C4 = C1/4;

fs = 2000;
dt = 1/fs;
t = 0:dt:6;

SD = 22;
MEAN = 220;
noise = normrnd(MEAN, SD, size(t));
P_in = noise; % Input P(t)
%%
y_save = zeros(6, length(t));
y = zeros(6, 1);
q_y = zeros(6,1);
%%
for n = 2:length(t)
    %%
    k = zeros(length(y), 4);
    for jr=1:4
        %%
        %%%% y(t + dt/2*k1), y(t + dt/2*k2), y(t + dt/2*k3)
        if jr ~= 1 
            y = y_save(:, n-1) + dt.*ra(jr).*k(:, jr-1);
        else
            y = y_save(:, n-1);
        end
        %%%% x(t + dt/2*k1), x(t + dt/2*k2), x(t + dt/2*k3)
        x_03 = sigmoidal_function(e0, R, v0, (y(3) - y(5)));
        x_14 = P_in(n-1) + C2.*sigmoidal_function(e0, R, v0, C1*y(1));
        x_25 = C4.*sigmoidal_function(e0, R, v0, C3*y(1));
        
        y_03 = y(1:2, 1);
        y_14 = y(3:4, 1);
        y_25 = y(5:6, 1);
        %% 4th-order Runge Kutta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        dy_03 = postsynaptic_potential_equation(y_03(2), y_03(1), A, a, x_03);
        dy_14 = postsynaptic_potential_equation(y_14(2), y_14(1), A, a, x_14);
        dy_25 = postsynaptic_potential_equation(y_25(2), y_25(1), B, b, x_25);
        
        k(:, jr) = [dy_03;dy_14;dy_25];
        r = dt/6.*(rb(jr).*k(:,jr));
        y = y + r;
    end
    %%
    y_save(:, n) = y;
end
y0 = y_save(1,:);
y3 = y_save(2,:);
y1 = y_save(3,:);
y4 = y_save(4,:);
y2 = y_save(5,:);
y5 = y_save(6,:);

eeg = y1-y2;
[f, fftp] = calc_fft_spectrum(eeg(t>=0.5), fs);

subplot(3,1,1)
plot(t, P_in)
xlabel('duration (s)')
ylabel('amplitude (a.u.)')
title('input signal')

subplot(3,1,2)
plot(t(2:end), eeg(2:end))
xlabel('duration (s)')
ylabel('amplitude (mV)')
title('simulated signal')

subplot(3,1,3)
plot(f, 10.*log(fftp));
xlim([0 50])
xlabel('frequency (Hz)')
ylabel('power (dB)')
title('power spectral density')