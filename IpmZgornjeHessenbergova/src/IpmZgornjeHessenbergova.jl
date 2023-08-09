module IpmZgornjeHessenbergova
using LinearAlgebra
import Base: \

export hessenberg

"""
    ZgornjiHessenberg{T}(H::Matrix{T})

Konstruktor za podatkovni tip `ZgornjiHessenberg` z zgornje Hessenbergovo matriko `H`.
"""
struct ZgornjiHessenberg{T <: Real}
    H::Matrix{T}
end


"""
    SpodnjaTridiagonalna{T}(data::Vector{T})

Konstruktor za podatkovni tip `SpodnjaTridiagonalna` z neničelnimi elementi pod diagonalo `data`.
"""
struct SpodnjaTridiagonalna{T <: Real}
    data::Vector{T}  # Neničelni elementi pod diagonalo
    function SpodnjaTridiagonalna{T}(data::Vector{T}) where T
        new{T}(data)
    end
end


"""
    hessenberg(A::Matrix{T}) -> ZgornjiHessenberg{T}, Matrix{T}

Izračuna zgornji Hessenbergov razcep matrike `A` in vrne zgornjo Hessenbergovo matriko ter matriko prehoda `Q`.
"""
function hessenberg(A::Matrix{T}) where T
    n = size(A, 1)
    Q = Matrix{T}(I, n, n)
    H = copy(A)

    for k in 1:n-2
        
        x = H[k+1:end, k]
        
        
        v = zeros(T, n-k)
        v[1] = sign(x[1])*norm(x, 2)
        
        v -= x
        
        
        if isequal(v, zeros(T, length(v)))
            continue
        end

        beta = 2 / norm(v)^2
        
        P = Matrix{T}(I, n-k, n-k) .- beta * (v * v')

        Q_k = Matrix{T}(I, n, n)
        Q_k[k+1:end, k+1:end] = P
        
        
        Q = Q_k * Q
        H = Q_k * H * Q_k'

    end

    return ZgornjiHessenberg(H), Q
end

"""
    lu(A::ZgornjiHessenberg{T}) -> SpodnjaTridiagonalna{T}, Matrix{T}, Matrix{Int}

Izvede LU razcep nad zgornjo Hessenbergovo matriko `A.H` in vrne razcepljene matrike `L`, `U` ter matriko permutacij `P`.
"""
function lu(A::ZgornjiHessenberg{T}) where T
    H = A.H
    n, m = size(H)
    v = zeros(T, n)
    U = copy(H)
    P = Matrix{Int}(I, n, n)  # Identiteta, ki sledi zamenjave vrstic

    for k in 1:n-1
        # Delno pivotarenje
        max_val = abs(U[k, k])
        max_index = k

        for i in k+1:n
            if abs(U[i, k]) > max_val
                max_val = abs(U[i, k])
                max_index = i
            end
        end

        if max_val < 1e-16
            error("Pivotni element je blizu nič. Potrebna uporaba druge metode.")
        end

        U[[k, max_index], k:end] = U[[max_index, k], k:end]
        P[[k, max_index], :] = P[[max_index, k], :]

        v[k+1] = U[k+1, k] / U[k, k]
        U[k+1, k:end] .-= v[k+1] * U[k, k:end]
    end

    v = v[end-(n-2):end]

    return SpodnjaTridiagonalna{T}(v), U, P
end

"""
    operator(L::SpodnjaTridiagonalna{T}, b::Vector{T}) where T

Operator za tip `SpodnjaTridiagonal`. Rešuje sistem enačb `L * y = x_n`.
"""
function \(L::SpodnjaTridiagonalna{T}, b::Vector{T}) where T
    n = length(L.data)
    x = copy(b)

    # Naprej
    for i in 1:n
        x[i+1] -= L.data[i] * x[i]
    end

    return x
end

"""
    inv_lastni(A::Matrix{T}, l::Real) -> Vector{T}, Real

Rešuje problem za prirejeno matriko polinoma, ki je podana z matriko `A` in začetnim približkom za lastno vrednost `l`.
"""
function inv_lastni(A::Matrix{T}, l::Real) where T
    zg_hess, sp_hess = hessenberg(A)
    H = zg_hess.H
    n = size(H, 1)
    lambda = l

    In = Matrix{Float64}(I, n, n)
    tol = 1e-16
    maxiter = 10^7
    
    # LU razcep matrike H - l*I
    L, U, P = lu(ZgornjiHessenberg(H - l*In))

    # Preverimo, ali je matrika U singularna
    if any(abs.(diag(U)) .< 1e-12)
        @info "Začetni približek za lastno vrednost je že lastna vrednost matrike."
        return lambda, NaN
    end

    # Nastavimo začetni lastni vektor
    xi = ones(Float64, n)
    xi = xi.*l

    for i in 1:maxiter
        # Rešimo sistem enačb L * y = x_n
        y = L \ xi
        
        # Rešimo sistem enačb U * x_{n+1} = y
        xi = U \ y

        # Normaliziramo xi
        xi /= norm(xi, 2)

        lambda = dot(xi, H*xi)

        # Preverimo konvergenco
        if norm((H - l*I)*xi, 2) < tol
            return xi, lambda
        end
    end

    # Metoda ni konvergirala
    return lambda, xi
end

end # module IpmZgornjeHessenbergova
