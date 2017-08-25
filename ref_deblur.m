function [ deblur, psf ] = ref_deblur( base, blur, psf_sz, cct_sz, alpha, beta, iter )
% Use the reference channel to restore blur channels
% Refer to ICCV'17 paper Revisiting Cross-channel Information Transfer ...
% for technical details

%% 0. Precompute 
sz = size(base);
deblur = zeros(sz(1), sz(2));
param_num = 4;
psf_step = ceil(psf_sz / 5);
cct_step = ceil(cct_sz / 5);

psf_num = ceil((sz(1) - psf_sz) / psf_step + 1) * ceil((sz(2) - psf_sz) / psf_step + 1);
cct_num = ceil((sz(1) - cct_sz) / cct_step + 1) * ceil((sz(2) - cct_sz) / cct_step + 1);

hamming_psf = repmat(reshape(hamming(psf_sz) * hamming(psf_sz)', psf_sz ^ 2, 1), [1, psf_num]);
hamming_psf_sum = col2imstep(hamming_psf, sz, [psf_sz, psf_sz], [psf_step, psf_step]);
hamming_cct = repmat(reshape(hamming(cct_sz) * hamming(cct_sz)', cct_sz ^ 2, 1), [1, cct_num]);
hamming_cct_sum = col2imstep(hamming_cct, sz, [cct_sz, cct_sz], [cct_step, cct_step]);

[X, Y] = meshgrid(1:param_num, 1:param_num);
i_index = reshape(bsxfun(@plus, repmat(Y, [1, 1, cct_num]), reshape(0:cct_num-1, 1 , 1, cct_num) * param_num), param_num^2, cct_num);
j_index = reshape(bsxfun(@plus, repmat(X, [1, 1, cct_num]), reshape(0:cct_num-1, 1 , 1, cct_num) * param_num), param_num^2, cct_num);

batch_base_psf_ft = fft2(reshape(im2colstep(base, [psf_sz, psf_sz], [psf_step, psf_step]), psf_sz, psf_sz, psf_num));
batch_blur_psf_ft = fft2(reshape(im2colstep(blur, [psf_sz, psf_sz], [psf_step, psf_step]), psf_sz, psf_sz, psf_num));
psf_priors = repmat(beta + alpha * (abs(fft2([1 -1], psf_sz, psf_sz)) .^ 2 + abs(fft2([1; -1], psf_sz, psf_sz)) .^ 2), [1, 1, psf_num]);

[base_gx, base_gy] = imgradientxy(base,'intermediate');

batch_base_mat = ones(cct_sz ^ 2, param_num, cct_num);
batch_base_mat(:, 2, :) = reshape(im2colstep(base, [cct_sz, cct_sz], [cct_step, cct_step]), cct_sz^2, 1, cct_num);
batch_base_mat(:, 3, :) = reshape(im2colstep(base_gx, [cct_sz, cct_sz], [cct_step, cct_step]), cct_sz^2, 1, cct_num);
batch_base_mat(:, 4, :) = reshape(im2colstep(base_gy, [cct_sz, cct_sz], [cct_step, cct_step]), cct_sz^2, 1, cct_num);

batch_blur = im2colstep(blur, [cct_sz, cct_sz], [cct_step, cct_step]);

%% 1. Initial guess
batch_deblur = bsxfun(@times, mean(batch_blur, 1) ./ mean(squeeze(batch_base_mat(:, 2, :)), 1), squeeze(batch_base_mat(:, 2, :)));
deblur = col2imstep(batch_deblur .* hamming_cct, sz, [cct_sz, cct_sz], [cct_step, cct_step]) ./ hamming_cct_sum;

for i = 1:iter
    %% 2. PSF estimation
    batch_deblur_psf_ft = fft2(reshape(im2colstep(deblur, [psf_sz, psf_sz], [psf_step, psf_step]), psf_sz, psf_sz, psf_num));
    batch_psf_ft = conj(batch_deblur_psf_ft) .* batch_blur_psf_ft ./ (conj(batch_deblur_psf_ft) .* batch_deblur_psf_ft + psf_priors);
    
    batch_blur_base = real(reshape(ifft2(batch_psf_ft .* batch_base_psf_ft), psf_sz^2, psf_num));
    blur_base = col2imstep(batch_blur_base .* hamming_psf, sz, [psf_sz, psf_sz], [psf_step, psf_step]) ./ hamming_psf_sum;
    
    %% 3. Cross channel transfer (CCT)
    [blur_base_gx, blur_base_gy] = imgradientxy(blur_base,'intermediate');
    
    batch_blur_base_mat = ones(cct_sz ^ 2, param_num, cct_num);
    batch_blur_base_mat(:, 2, :) = reshape(im2colstep(blur_base, [cct_sz, cct_sz], [cct_step, cct_step]), cct_sz^2, 1, cct_num);
    batch_blur_base_mat(:, 3, :) = reshape(im2colstep(blur_base_gx, [cct_sz, cct_sz], [cct_step, cct_step]), cct_sz^2, 1, cct_num);
    batch_blur_base_mat(:, 4, :) = reshape(im2colstep(blur_base_gy, [cct_sz, cct_sz], [cct_step, cct_step]), cct_sz^2, 1, cct_num);
    
    A = reshape(mmx('square', batch_blur_base_mat, [], 't'), param_num^2, cct_num);
    b = mmx('mult', batch_blur_base_mat, reshape(batch_blur, cct_sz ^ 2, 1, cct_num), 'tn');
    AS = sparse(i_index, j_index, A, param_num * cct_num, param_num * cct_num);
    batch_coefficient = reshape(AS \ reshape(b, param_num * cct_num, 1), param_num, 1, cct_num);
    batch_deblur = squeeze(mmx('mult', batch_base_mat, batch_coefficient));
    deblur = col2imstep(batch_deblur .* hamming_cct, sz, [cct_sz, cct_sz], [cct_step, cct_step]) ./ hamming_cct_sum;
    
end
psf_half = 24;
kernel = real(ifft2(batch_psf_ft));
psf = kernel([end-psf_half+1:end, 1:psf_half+1], [end-psf_half+1:end, 1:psf_half+1], :);
end

