function [avg_err] = KLT(img, survive, bc_sz)
%KLT function that reads the image and reconstruct the image using KLT
%   INPUT
%   img: the system path of Image (recommand to use a square image)
%        e.g. "./lena.png"
%   survive: how many components to survive
%   bc_sz: the partitioned block size
%               e.g. if want to use 4*4 block, set block_size = 4
%   OUTPUT
%   avg_err: the squared error with current configuration

    im = imread(img); % read the image
    im = rgb2gray(im); % trans to grayscale
    imshow(im); % show image
    sz = size(im); % the image size (for lena, 512, 512)
    if mod(sz(1), bc_sz) % the block size need to be exact divided by full size
        error('ERROR. Input parameter should be divided by %.f', sz(1));
    end
    pcs = sz(1)/bc_sz; 
    pcs_pos = ones(1, pcs)*bc_sz;
    part_im = mat2cell(im, pcs_pos, pcs_pos); % partition the image
    bc_num = pcs^2; % total number of sub-images
    vector_len = bc_sz^2; % the vector length

    cat_mat = zeros([vector_len, bc_num]); % cat_mat is the concated matrix for average calc
    tmp = 1;
    for i = 1:pcs
        for j = 1:pcs
            tmp1 = reshape(transpose(part_im{i,j}), [vector_len 1]); % concat the subimg
            cat_mat(:,tmp) = tmp1;
            tmp = tmp + 1;
        end
    end
    x_avg = mean(cat_mat.'); % the average vector

    acc_cor = zeros([vector_len, vector_len]); % to save the accumulated corelevance matrix
    for i=1:bc_num
        tmp1 = cat_mat(:,i) - x_avg.'; % (vector)-(average vector)
        tmp2 = tmp1*transpose(tmp1);
        acc_cor = acc_cor + tmp2;
    end
    cor = 1/(bc_num - 1)*acc_cor; % the corelevance matrix

    [V,D] = eig(cor); % calculate the eigen vectors and eigenvalue.
    [~, ind_eig] = sort(diag(D)); % sort, from low to high
    ind_eig = ind_eig(vector_len-survive+1:end); % get index of the sorted eigen value, for further drop

    bc_err = zeros([bc_num, 1]);
    for j = 1:bc_num
        curr = cat_mat(:, j);
        est = zeros(vector_len, 1); % to save the estimated vector
        for i = 1:survive % drop some value till survive number survived
            ind = ind_eig(i); % the eigen value's index
            est = est + transpose(curr)*V(:, ind)*V(:, ind); % calculate the est for each direction
        end
        err = (est - curr).^2; % squared err
        bc_err(j) = sum(err);
    end
    avg_err = mean(bc_err); % average error for all blocks
end
