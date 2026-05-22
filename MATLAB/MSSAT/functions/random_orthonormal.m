function I = random_orthonormal(n, m)

    if n < m
        warning('dimension is lower than the number of vectors');
        I = [];
        return;
    end

    I = zeros(n, m);

    I(:, 1) = randn(n, 1);
    I(:, 1) = I(:, 1) / norm(I(:, 1));

    for i = 2:m
        vec = randn(n, 1);

        proj = I(:, 1:(i-1)) * (I(:, 1:(i-1))' * vec);
        vec = vec - proj;

        I(:, i) = vec / norm(vec);
    end
end
