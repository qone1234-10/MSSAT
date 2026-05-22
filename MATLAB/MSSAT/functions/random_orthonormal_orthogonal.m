function I = random_orthonormal_orthogonal(m, mat)

    n = size(mat, 1);

    if n < m
        warning('dimension is lower than the number of vectors');
        I = [];
        return;
    end

    I = zeros(n, m);

    vec = randn(n, 1);
    I(:, 1) = vec / norm(vec);

    I = I - mat / (mat' * mat) * mat' * I;

    I(:, 1) = I(:, 1) / norm(I(:, 1));

    for i = 2:m
        vec = randn(n, 1);
        vec = vec -  mat / (mat' * mat) * mat' * vec;

        proj = I(:, 1:(i-1)) * I(:, 1:(i-1))' * vec;
        vec = vec - proj;

        I(:, i) = vec / norm(vec);
    end
end
