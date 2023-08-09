using LinearAlgebra, IpmZgornjeHessenbergova

A = Float64[1.0 0.0 2.0 3.0;
            -1.0 0.0 5.0 2.0;
            0.0 -2.0 0.0 0.0;
            0.0 0.0 2.0 0.0]

# Primer uporabe metode hessenberg,
# ki naredi hessenbergov razcep
# zg_hess je tipa ZgornjiHessenberg, marika Q je tipa Matrix
zg_hess, Q = IpmZgornjeHessenbergova.hessenberg(A)

# Zgornje hessenbergova matrika 
H = zg_hess.H

#Primer LU razcepa zgornje Hessenbergove matrike 
#matrika L je spodnjetrikotna, matrika U je zgornjetrikotna in P je projekcijska matrika
L, U, P = IpmZgornjeHessenbergova.lu(IpmZgornjeHessenbergova.ZgornjiHessenberg(A))

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

# Koeficienti polinom x^2 + 4x + 4
coeffs = [4.0, 4.0, 1.0]

# Izračun matrike polinoma
println("Prirejena matrika polinomu x^2 + 4x + 4:")
B = priredi_matriko_polinomu(coeffs)
println(B)
println()

# Izračun lastnih vrednosti z metodo eigen
println("Lastne vrednosti prirejene matrike polinoma x^2 + 4x + 4 z metodo eigen:")
eigenvals = eigen(B).values
println(eigenvals)
println()

# Izračun lastnih vrednosti z metodo inv_lastni
println("Lastne vrednosti prirejene matrike polinoma x^2 + 4x + 4 z metodo inv_lastni:")
lastna_vrednost, lastni_vektor = IpmZgornjeHessenbergova.inv_lastni(B, 2)
println(lastna_vrednost)
println()

# Koeficienti polinom (x-1)*(x-2)*(x-3) = x^3 - 6x^2 + 11x - 6
coeffs = [-6.0, 11.0, -6.0, 1.0]

# Izračun matrike polinoma
println("Prirejena matrika polinomu x^3 - 6x^2 + 11x - 6:")
B = priredi_matriko_polinomu(coeffs)
println(B)
println()

# Izračun lastnih vrednosti z metodo eigen
println("Lastne vrednosti prirejene matrike polinoma x^3 - 6x^2 + 11x - 6 z metodo eigen:")
eigenvals = eigen(B).values
println(eigenvals)
println()

# Izračun lastnih vrednosti z metodo inv_lastni z začtenim približkom za lastno vrednost 2.0
println("Lastne vrednosti prirejene matrike polinoma x^3 - 6x^2 + 11x - 6 z metodo inv_lastni za lastno vrednost 2.0:")
lastna_vrednost, lastni_vektor = IpmZgornjeHessenbergova.inv_lastni(B, 2.0)
println(lastna_vrednost)
println()

# Izračun lastnih vrednosti z metodo inv_lastni z začtenim približkom za lastno vrednost 6.0
println("Lastne vrednosti prirejene matrike polinoma x^3 - 6x^2 + 11x - 6 z metodo inv_lastni za lastno vrednost 6.0:")
lastna_vrednost, lastni_vektor = IpmZgornjeHessenbergova.inv_lastni(B, 6.0)
println(lastna_vrednost)
println()