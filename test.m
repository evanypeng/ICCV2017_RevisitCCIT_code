%% *** For research purposes only ***

% A simplified, non-accelerated demo implementation of the algorithm described in: 
% Revisiting Cross-channel Information Transfer for Chromatic Aberration Correction
% Accepted to ICCV 2017 (IEEE)
% Authors: Tiancheng Sun, Yifan Peng, Wolfgang Heidrich
% If you found any bugs in the demo implementation, please contact us at 
% kevin.kingo0627@gmail.com / evanpeng@cs.ubc.ca  
% We appreciate that.

%% Parameters setting for aberration correction
psf_sz = 500;  % PSF estimation window size
win_sz = 60;   % CCT implementation window size
alpha = 0.3;
beta = 0.3;
iter = 3;

use_blind = 0;

%% Run the demo experiment
path = 'im/';

im = im2double(imread([path, 'blur.png']));
deblurred = zeros(size(im));
deblurred(:, :, 2) = im(:, :, 2);

if use_blind
    filename = [path, 'green_blinddeconv.png'];
    if exist(filename, 'file')
        deblurred(:, :, 2) = im2double(imread(filename));
    end
end

[deblurred(:, :, 1), psf_red] = ref_deblur(deblurred(:, :, 2), im(:, :, 1), psf_sz, win_sz, alpha, beta, iter);
[deblurred(:, :, 3), psf_blue] = ref_deblur(deblurred(:, :, 2), im(:, :, 3), psf_sz, win_sz, alpha, beta, iter);

imwrite(deblurred, [path, 'ours_', num2str(psf_sz), '_', num2str(win_sz), '.jpg']);

