%% Parameters
highRes = 2048;         % High resolution image size (assumed square)
movieSize = 256;        % Size of each movie frame (sampling window)
numFrames = 100;        % Number of movie frames
D = 1;                  % Diffusion constant (pixels^2 per frame)

%% Generate a high-resolution image with 1/f amplitude scaling
% Create frequency grid (with zero-centered frequencies)
[u, v] = meshgrid(-highRes/2:highRes/2-1, -highRes/2:highRes/2-1);
rho = sqrt(u.^2 + v.^2);
% Avoid division by zero at the DC component
rho(rho == 0) = 1;  
% For a power spectrum ~1/f^2, the amplitude spectrum scales as 1/f.
amplitude = 1 ./ rho;

% Create a random phase uniformly distributed between 0 and 2pi
phi = 2*pi*rand(highRes, highRes);

% Construct the Fourier domain representation
FT = amplitude .* exp(1i*phi);
% Optionally, set the DC component to zero (remove overall bias)
% FT(highRes/2+1, highRes/2+1) = 0;

% Inverse Fourier transform to get the spatial-domain image.
% Use ifftshift to put the zero-frequency component back in place.
imgHighRes = real(ifft2(ifftshift(FT)));

% Compute the Fourier transform and power spectrum of the static image
F_static = fftshift(fft2(imgHighRes));
P_static = abs(F_static).^2;

% Create spatial frequency grid (in cycles per pixel)
[u, v] = meshgrid(-highRes/2:highRes/2-1, -highRes/2:highRes/2-1);
rho = sqrt(u.^2 + v.^2);

% Radial binning of the power spectrum
numBins = 50;
r_max = max(rho(:));
binEdges = linspace(0, r_max, numBins+1);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
P_radial = zeros(1, numBins);

% Average the power over annular bins.
for i = 1:numBins
    mask = (rho >= binEdges(i)) & (rho < binEdges(i+1));
    if any(mask(:))
        P_radial(i) = mean(P_static(mask));
    else
        P_radial(i) = NaN;
    end
end

%% Theoretical power spectrum: natural images are expected to have a power law ~1/k^2
% Create the theoretical curve using the bin center frequencies.
P_theo = 1 ./ (binCenters.^2);
% Scale the theoretical curve to the empirical data.
% (Since the absolute scale is arbitrary, we match the maximum values.)
% scale_factor = max(P_radial(~isnan(P_radial))) / max(P_theo(~isnan(P_theo)));
% P_theo = scale_factor * P_theo;

% Plot the results on a log-log plot for clarity
figure;
loglog(binCenters, P_radial, 'b-', 'LineWidth', 2);
hold on;
loglog(binCenters, P_theo, 'r--', 'LineWidth', 2);
xlabel('Spatial Frequency (cycles/pixel)');
ylabel('Power');
legend('Empirical', 'Theoretical (1/k^2)');
title('Radial Average of Static Image Power Spectrum vs. Theoretical 1/k^2');
grid on;


%% Precompute grid for interpolation on the high-res image
[xHigh, yHigh] = meshgrid(1:highRes, 1:highRes);

% Simulate random walk for subpixel motion
% Start from the center of the high-res image
xPos = zeros(numFrames,1);
yPos = zeros(numFrames,1);
xPos(1) = highRes/2;
yPos(1) = highRes/2;

% For each subsequent frame, update the position with a Gaussian step.
for t = 2:numFrames
    dx = sqrt(2*D) * randn;
    dy = sqrt(2*D) * randn;
    xPos(t) = xPos(t-1) + dx;
    yPos(t) = yPos(t-1) + dy;
end

% Prepare grid for the movie frame sampling
% Define the local coordinate grid for the movie window.
halfWin = movieSize/2;
% Create a grid so that the window is centered; using subpixel offsets.
[xLocal, yLocal] = meshgrid(linspace(-halfWin+0.5, halfWin-0.5, movieSize), ...
                            linspace(-halfWin+0.5, halfWin-0.5, movieSize));

% Create movie frames by sampling the high-res image using interp2
frames(numFrames) = struct('rawdata', [],'cdata', [], 'colormap', []);
for t = 1:numFrames
    % Calculate the coordinates in the high-res image where the movie window is sampled.
    Xq = xPos(t) + xLocal;
    Yq = yPos(t) + yLocal;
    % Use linear interpolation; pixels falling outside the image are set to 0.
    frame = interp2(xHigh, yHigh, imgHighRes, Xq, Yq, 'linear', 0);
    % Convert to an 8-bit image for display (adjust as desired)
    frameReadyForInt8 = frame/max(frame(:));
    frames(t).rawdata = frame;
    frames(t).cdata = im2uint8(frameReadyForInt8);
    frames(t).colormap = colormap('default');
end

%% Play the movie
figure;
movie(gcf, frames, 1, 10);  % Play once at 10 frames per second

%% Optionally, write the movie to a video file
% Uncomment the following section to save the movie as an .avi file.
%{
v = VideoWriter('subpixelRandomWalk.avi');
v.FrameRate = 10;
open(v);
for t = 1:numFrames
    writeVideo(v, frames(t).cdata);
end
close(v);
%}




%% Parameters & Data Setup
% (Assumes 'frames' from the previous movie-generation code is in your workspace)
movieSize = size(frames(1).rawdata, 1);  % e.g., 256
numFrames = length(frames);            % e.g., 100

% Convert movie frames (stored as 8-bit images) into a 3D double array
movieData = zeros(movieSize, movieSize, numFrames);
for t = 1:numFrames
    % Normalize to [0,1]
    movieData(:,:,t) = double(frames(t).rawdata);
end

%% Compute 3D FFT and Empirical Spatiotemporal Power Spectrum
% Normalize by the total number of samples to get the correct amplitude
N_total = movieSize * movieSize * numFrames;
F = fftshift(fftn(movieData)) / N_total;
P_emp = abs(F).^2;

%% Frequency Axes
% Spatial frequency axis (cycles per pixel). We assume pixel spacing = 1.
fx = (-movieSize/2:movieSize/2-1) / movieSize;
fy = fx;
[FX, FY] = meshgrid(fx, fy);
fr = sqrt(FX.^2 + FY.^2);  % Radial spatial frequency (cycles per pixel)

% Temporal frequency axis (cycles per frame)
ft = (-numFrames/2:numFrames/2-1) ;


% Create 2D spatial frequency grid for one time slice.
[FX_2, FY_2] = meshgrid(fx, fy);
% Replicate the spatial grids along the temporal dimension.
FX_3 = repmat(FX_2, [1, 1, numFrames]);
FY_3 = repmat(FY_2, [1, 1, numFrames]);

% Replicate the temporal frequency axis across the spatial dimensions.
FT_3 = reshape(ft, [1, 1, numFrames]);
FT_3 = repmat(FT_3, [movieSize, movieSize, 1]);

% Compute the radial spatial frequency for each point in the 3D grid.
F_rad = sqrt(FX_3.^2 + FY_3.^2);

% Empirical Temporal Spectrum for Spatial Frequency Bins
% Average the spatiotemporal power over annular bins in spatial frequency.
% Define spatial frequency bin edges (in cycles/pixel)
binEdges = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5];
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
numBins = length(binCenters);

% Preallocate: rows = spatial frequency bin, columns = temporal frequency samples
P_emp_temp = zeros(numBins, numFrames);

% For each bin, average the power (over spatial indices) for every temporal frequency slice.
for bt = 1:numBins
    % Create a mask selecting spatial frequencies in the current bin.
    mask = (fr >= binEdges(bt)) & (fr < binEdges(bt+1));
    for it = 1:numFrames
        slice = P_emp(:,:,it);
        if any(mask(:))
            P_emp_temp(bt, it) = mean(slice(mask));
        end
    end
end

%% --- 3D Theoretical Spectrum ---

% Compute the theoretical 3D spatiotemporal power spectrum.
% Using the formula:
%   P_theo(f, ft) = 2D / [ (2π f_t)^2 + (D*(2π f)^2)^2 ]
P_theo = 2*D ./ ( (2*pi*FT_3).^2 + (D*(2*pi*F_rad).^2).^2 );

% --- Binning the Theoretical Spectrum ---

% Preallocate binned theoretical temporal spectra.
P_theo_temp = zeros(numBins, numFrames);

% For each annular spatial frequency bin, average the theoretical power over the
% spatial dimensions for each temporal frequency slice.
for bt = 1:numBins
    % Create a mask selecting spatial frequencies in the current bin.
    mask = (fr >= binEdges(bt)) & (fr < binEdges(bt+1));
    for it = 1:numFrames
        slice = P_theo(:,:,it);
        if any(mask(:))
            P_theo_temp(bt, it) = mean(slice(mask));
        end
    end
end

%% --- (Optional) Rescale Theoretical Spectrum ---
% As the theoretical derivation is only defined up to a multiplicative constant,
% you might want to match one of the bins (for example, the first) to the empirical data.
 % alpha = max(P_emp_temp(1,:)) / max(P_theo_temp(1,:));
 % P_theo_temp = alpha * P_theo_temp;

%% --- Plotting: Compare Empirical and 3D Theoretical Temporal Spectra ---

figure;
for bt = 1:numBins
    subplot(numBins, 1, bt);
    plot(ft, P_emp_temp(bt, :)./sum(P_emp_temp(bt, :)), 'b', 'LineWidth', 1.5); hold on;
    plot(ft, P_theo_temp(bt, :)./sum(P_theo_temp(bt, :)), 'r--', 'LineWidth', 1.5);
    xlabel('Temporal Frequency (cycles/frame)');
    ylabel('Power');
    title(sprintf('Spatial Frequency Bin: [%.3f, %.3f] (Center = %.3f)', ...
          binEdges(bt), binEdges(bt+1), binCenters(bt)));
    legend('Empirical', 'Theoretical');
end
sgtitle('Comparison of Empirical vs. 3D Theoretical Temporal Power Spectra');