% Grid parameters
N = 256;   % Number of points in each dimension
x = linspace(-5, 5, N);
y = linspace(-5, 5, N);
[X, Y] = meshgrid(x, y);

% Example of a 2D even, real-valued signal (Gaussian)
sigma = .1;  % Standard deviation
f = exp(-(X.^2 + Y.^2)/(2*sigma^2));

% Compute the 2D FFT of the signal
F = fft2(f);
F1=F;
% Ensure the FFT is real (using fftshift to center the zero frequency component)
F_shifted = fftshift(F);
F_real = real(F_shifted);
F_imag = imag(F_shifted);

% Plot the original signal
figure;
subplot(3, 3, 1);
imagesc(x, y, f);
title('Original 2D Signal (Gaussian)');
xlabel('x');
ylabel('y');
colorbar;

% Plot the real part of the Fourier transform
subplot(3, 3, 2);
imagesc(x, y, F_real);
title('Real Part of 2D FFT');
xlabel('kx');
ylabel('ky');
colorbar;

% Plot the imaginary part of the Fourier transform
subplot(3, 3, 3);
imagesc(x, y, F_imag);
title('Imaginary Part of 2D FFT');
xlabel('kx');
ylabel('ky');
colorbar;


% Check if imaginary part is near zero (numerical stability)
disp('Maximum absolute value of imaginary part of FFT:');
disp(max(max(abs(imag(F_shifted)))));



N = 256;   % Number of points in each dimension
x = (2*(0:N-1)/N - 1)*5;
y = (2*(0:N-1)/N - 1)*5;
[X, Y] = meshgrid(x, y);

% Example of a 2D even, real-valued signal (Gaussian)
sigma = .1;  % Standard deviation
f = exp(-(X.^2 + Y.^2)/(2*sigma^2));

% Compute the 2D FFT of the signal
F = fft2(f);
F2=F;
% Ensure the FFT is real (using fftshift to center the zero frequency component)
F_shifted = fftshift(F);
F_real = real(F_shifted);
F_imag = imag(F_shifted);

% Plot the original signal

subplot(3, 3, 3+1);
imagesc(x, y, f);
title('Original 2D Signal (Gaussian)');
xlabel('x');
ylabel('y');
colorbar;

% Plot the real part of the Fourier transform
subplot(3, 3, 3+2);
imagesc(x, y, F_real);
title('Real Part of 2D FFT');
xlabel('kx');
ylabel('ky');
colorbar;

% Plot the imaginary part of the Fourier transform
subplot(3, 3, 3+3);
imagesc(x, y, F_imag);
title('Imaginary Part of 2D FFT');
xlabel('kx');
ylabel('ky');
colorbar;


% Check if imaginary part is near zero (numerical stability)
disp('Maximum absolute value of imaginary part of FFT:');
disp(max(max(abs(imag(F_shifted)))));




% Check if imaginary part is near zero (numerical stability)
disp('Maximum absolute value of imaginary part of FFT:');
disp(max(max(abs(imag(F_shifted)))));



N = 256;   % Number of points in each dimension
x = (2*(0:N-1)/N - 1)*5;
y = (2*(0:N-1)/N - 1)*5;
[X, Y] = meshgrid(x, y);

% Example of a 2D even, real-valued signal (Gaussian)
sigma = .1;  % Standard deviation
f = exp(-(X.^2 + Y.^2)/(2*sigma^2));
f = fftshift(f);

% Compute the 2D FFT of the signal
F = fft2(f);
F2=F;
% Ensure the FFT is real (using fftshift to center the zero frequency component)
F_shifted = fftshift(F);
F_real = real(F_shifted);
F_imag = imag(F_shifted);

% Plot the original signal

subplot(3, 3, 6+1);
imagesc(x, y, f);
title('Original 2D Signal (Gaussian)');
xlabel('x');
ylabel('y');
colorbar;

% Plot the real part of the Fourier transform
subplot(3, 3, 6+2);
imagesc(x, y, F_real);
title('Real Part of 2D FFT');
xlabel('kx');
ylabel('ky');
colorbar;

% Plot the imaginary part of the Fourier transform
subplot(3, 3, 6+3);
imagesc(x, y, F_imag);
title('Imaginary Part of 2D FFT');
xlabel('kx');
ylabel('ky');
colorbar;


% Check if imaginary part is near zero (numerical stability)
disp('Maximum absolute value of imaginary part of FFT:');
disp(max(max(abs(imag(F_shifted)))));