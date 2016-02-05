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

function square(N)
    l = [((-1, 0), 1), ((0, -1), 1), ((0, 0), -4), ((0, 1), 1), ((1, 0), 1)]
end

function hexagonal(N)
    l = [((-1, 0), 1), ((0, -1), 1), ((0, 0), -6), ((0, 1), 1), ((1, 0), 1), ((1, 1), 1), ((-1, -1), 1)]
end

function boundary(ij, lower, upper)
    i, j = ij
    (i < lower || upper < i || j < lower || upper < j)
end

function k2ij(k, N)
    k = k - 1
    floor(Int, k/N)+1, k%N+1
end

function ij2k(ij, N)
    i, j = ij
    N*(i-1) + j
end

# 1 based
# u_{i,j} = u_{1,1}, u_{1,2}, ...
# k = 1, 2, ...
function make_A(N, f)
    A = zeros(N^2, N^2)
    for k in 1:N^2
#        print(k, ": ")
#        print(k2ij(k, N))
#        print("\n")
        i, j = k2ij(k, N)
        nn_ijv_list = map(ijv -> (map(+, (i,j), ijv[1]), ijv[2]), f(N))
        nn_ijv_list = filter(ijv -> !boundary(ijv[1], 1, N), nn_ijv_list)
        nn_k_list = map(ijv -> (ij2k(ijv[1], N), ijv[2]), nn_ijv_list)
        for (kk,v) in nn_k_list
            A[k, kk] = v
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
