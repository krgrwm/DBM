#A = Float64[
#-4   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0;
# 1  -4   1   0   0   1   0   0   0   0   0   0   0   0   0   0;
# 0   1  -4   1   0   0   1   0   0   0   0   0   0   0   0   0;
# 0   0   1  -4   0   0   0   1   0   0   0   0   0   0   0   0;
# 1   0   0   0  -4   1   0   0   1   0   0   0   0   0   0   0;
# 0   1   0   0   1  -4   1   0   0   1   0   0   0   0   0   0;
# 0   0   1   0   0   1  -4   1   0   0   1   0   0   0   0   0;
# 0   0   0   1   0   0   1  -4   0   0   0   1   0   0   0   0;
# 0   0   0   0   1   0   0   0  -4   1   0   0   1   0   0   0;
# 0   0   0   0   0   1   0   0   1  -4   1   0   0   1   0   0;
# 0   0   0   0   0   0   1   0   0   1  -4   1   0   0   1   0;
# 0   0   0   0   0   0   0   1   0   0   1  -4   0   0   0   1;
# 0   0   0   0   0   0   0   0   1   0   0   0  -4   1   0   0;
# 0   0   0   0   0   0   0   0   0   1   0   0   1  -4   1   0;
# 0   0   0   0   0   0   0   0   0   0   1   0   0   1  -4   1;
# 0   0   0   0   0   0   0   0   0   0   0   1   0   0   1  -4
# ]
# このAは間違っている

function out_of_bounds(i, j, N)
    (i <= 0 || N < i || j <= 0 || N < j)
end

function square(N)
    l = [((-1, 0), 1), ((0, -1), 1), ((0, 0), -4), ((0, 1), 1), ((1, 0), 1)]
end

function hexagonal(N)
    l = [((-1, 0), 1), ((0, -1), 1), ((0, 0), -6), ((0, 1), 1), ((1, 0), 1), ((1, 1), 1), ((-1, -1), 1)]
end

# 0が基準
function displacement(i, j, N)
    N*i + j
end

function disp_list(N, f)
    map(x->(displacement(x[1]..., N), x[2]), f(N))
end
disp_list(5, hexagonal)

function make_A(N, f)
    A = zeros(N^2, N^2)
    dl = disp_list(N, f)
    for i in 1:N^2
        for (d,v) in dl
            if !out_of_bounds(i, i+d, N^2)
                A[i, i+d] = v
            end
        end
    end
    A
end
    
A = make_A(10, square)
A = make_A(10, hexagonal)

function check_convergence(M)
    eigvals_check = !(false in map(x->abs(x) < 1.0, eigvals(M)))
    Dict(
    "1"           => norm(M, 1)   < 1.0,
    "2"           => norm(M, 2)   < 1.0,
    "Inf"         => norm(M, Inf) < 1.0,
    "Frobenius"   => vecnorm(M)   < 1.0,
    "eigenvalues" => eigvals_check
    )
end

U = triu(A, 1)
L = tril(A, -1)
D = A |> diag |> diagm
# jacobi
M_jacobi = -inv(D) * (L + U)
check_convergence(M_jacobi)

# gauss-sidel
M_gauss_seidel = -inv(D+L) * U
check_convergence(M_gauss_seidel)

# SOR
omega = 1.5
M_sor = inv(D+omega*L) * ((1-omega)*D - omega*U)
check_convergence(M_sor)
