<img width="1919" height="1077" alt="image" src="https://github.com/user-attachments/assets/d2978fd9-a547-4666-8e08-b2177b7affa8" /># Design-of-FIR-Filters-using-hanning-window
REG N0:212223060247
#DESIGN OF FIR DIGITAL FILTER 

# AIM: 
          
  To generate design of low pass FIR digital filter using SCILAB 

# APPARATUS REQUIRED: 

  PC Installed with SCILAB 

# PROGRAM 
LPF
```
// FIR Low-Pass Filter using Fourier Series
// Hanning window, cutoff = 0.5π radians/sample, N = 21

clear;
clc;

// Filter specs
fs = 1;              // Normalized sampling frequency
fc = 0.25;           // Normalized cutoff (cycles/sample)
N = 21;              // Filter length (odd)

// Time index centered around zero
n = -(N-1)/2 : (N-1)/2;

// Ideal impulse response (sinc)
h = zeros(1, N);
for i = 1:N
    if n(i) == 0 then
        h(i) = 2 * fc;
    else
        h(i) = sin(2 * %pi * fc * n(i)) / (%pi * n(i));
    end
end

// Apply Hanning window
w = 0.5 - 0.5 * cos(2 * %pi * (0:N-1) / (N-1));
h_win = h .* w;

// Display coefficients
disp('FIR LPF Coefficients (Hanning Window):');
disp(h_win);

// Frequency response
scf(0);
clf();

[Hf, fr] = frmag(h_win, 512);
subplot(2,1,1);
plot(fr * fs, Hf);
xtitle('Magnitude Response (Hanning Window)', 'Frequency (Hz)', 'Magnitude');

// Phase response
w = 2 * %pi * fr;
H_sym = poly(h_win, 'z');
H_eval = horner(H_sym, exp(-%i * w));
phase_rad = atan(imag(H_eval) ./ real(H_eval));
phase_deg = phase_rad * 180 / %pi;

subplot(2,1,2);
plot(fr * fs, phase_deg);
xtitle('Phase Response (Hanning Window)', 'Frequency (Hz)', 'Phase (degrees)');

// Filter info
disp('FIR LPF Design Complete');
disp('Cutoff Frequency: 0.5π radians/sample');
disp('Filter Length: ' + string(N));
```
HPF
```
clear;
clc;

// Filter specs
fs = 1;              // Normalized sampling frequency
fc = 0.25;           // Normalized cutoff (cycles/sample)
N = 21;              // Filter length (odd)

// Time index centered around zero
n = -(N-1)/2 : (N-1)/2;

// Ideal LPF impulse response (sinc)
h_lp = zeros(1, N);
for i = 1:N
    if n(i) == 0 then
        h_lp(i) = 2 * fc;
    else
        h_lp(i) = sin(2 * %pi * fc * n(i)) / (%pi * n(i));
    end
end

// Spectral inversion to get HPF
h_hp = -h_lp;
h_hp((N+1)/2) = 1 - h_lp((N+1)/2);  // center tap adjustment

// Apply Hanning window
w = 0.5 - 0.5 * cos(2 * %pi * (0:N-1) / (N-1));
h_win = h_hp .* w;

// Display coefficients
disp('FIR HPF Coefficients (Hanning Window):');
disp(h_win);

// Frequency response
scf(0);
clf();

[Hf, fr] = frmag(h_win, 512);
subplot(2,1,1);
plot(fr * fs, Hf);
xtitle('Magnitude Response (Hanning Window)', 'Frequency (Hz)', 'Magnitude');

// Phase response
w = 2 * %pi * fr;
H_sym = poly(h_win, 'z');
H_eval = horner(H_sym, exp(-%i * w));
phase_rad = atan(imag(H_eval) ./ real(H_eval));
phase_deg = phase_rad * 180 / %pi;

subplot(2,1,2);
plot(fr * fs, phase_deg);
xtitle('Phase Response (Hanning Window)', 'Frequency (Hz)', 'Phase (degrees)');

// Filter info
disp('FIR HPF Design Complete');
disp('Cutoff Frequency: 0.5π radians/sample');
disp('Filter Length: ' + string(N));
```

# OUTPUT
LPF
<img width="1221" height="480" alt="image" src="https://github.com/user-attachments/assets/1ea3989f-40f7-4146-93b2-d0c433b4c284" />
<img width="1221" height="480" alt="image" src="https://github.com/user-attachments/assets/d2f993f9-10f0-45e4-a9f9-5e90fc9be9ea" />

HPF
<img width="1238" height="418" alt="image" src="https://github.com/user-attachments/assets/6817b5a5-b6a2-4d9b-928e-af20e2678ca5" />

<img width="1918" height="1078" alt="image" src="https://github.com/user-attachments/assets/d1a0b13e-71af-47f3-9b3d-d75e853404eb" />

# RESULT
Thus FIR Low pass filter  High pass filter design using Hanning window is performed.
