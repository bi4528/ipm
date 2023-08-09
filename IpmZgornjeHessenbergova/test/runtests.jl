using Test, IpmZgornjeHessenbergova
using LinearAlgebra
using Polynomials

tolerance = 1e-10
n = 5
A = rand(n ,n)
zg_hess, Q = IpmZgornjeHessenbergova.hessenberg(A)
H = zg_hess.H
L, U, P = IpmZgornjeHessenbergova.lu(zg_hess)

#test če so ta pravi elementi matrike H ničleni elementi
function je_zgornje_hessenbergova(H::Matrix{T}, tol::Float64) where T
    n = size(H, 1)
    for i in 1:n, j in 1:n
        if j < i - 2 && abs(H[i, j]) > tol
            return false
        end
    end
    return true
end
@test je_zgornje_hessenbergova(H, tolerance)

#test če je matrika H zgornje Hessenbergova
@test all(abs.(Q*A*Q' - H) .< tolerance)

#test metode inv_lastni
function ustvari_polinom_iz_korenov(roots::Vector{T}) where T
    poly = fromroots(roots)  # Create a polynomial with specified roots
    return poly.coeffs
end
roots = rand(n)
roots = round.(roots .* 10.0) .+ 1.0
roots[end] = 1.0
coeffs = ustvari_polinom_iz_korenov(roots)

# Funkcija za izračun matrike polinoma
function priredi_matriko_polinomu(coeffs::Vector{Float64})
    n = length(coeffs) - 1
    A = zeros(n, n)
    for i in 1:n-1
        A[i+1, i] = 1.0
    end
    A[:, end] = -coeffs[1:end-1]
    
    return A
end

B = priredi_matriko_polinomu(coeffs)

lastna_vrednost, lastni_vektor = IpmZgornjeHessenbergova.inv_lastni(B, -1.0)

@test any(abs.(roots .- lastna_vrednost) .< tolerance)




