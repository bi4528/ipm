# Modul `IpmZgornjeHessenbergova`

Modul `IpmZgornjeHessenbergova` vsebuje implementacijo inverzne potenčne metode z uporabo zgornje Hessenbergove matrike, ki jo lahko uporabimo za iskanje ničle polinoma:
$x^n + a_{n-1}x^{n-2}+...+a_1x+a_0$  
Vsakemu takšnemu polinomu lahko priredimo matriko 
```math
\begin{bmatrix}
0 & 0 & \cdots & 0 & -a_0 \\
1 & 0 & \cdots & 0 & -a_1 \\
0 & 1 & \cdots & 0 & -a_2 \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & \cdots & 1 & -a_{n-1}
\end{bmatrix}
```
in z inverzno potenčno metodo lahko pridobimo njeno lastno vrednost, ki je ničla tega polinoma. 
Inverzna potenčna metoda je implementirana tako, da se na poljubni matriki z metodo `hessenberg` najprej naredi Hessenbergov razcep iz katerega pridobimo zgornje Hessenbergovo matriko. Potem pa naredimo LU razcep pridobljene zgornje Hessenbergove matrike z uporabo metode `lu` iz modula `IpmZgornjeHessenbergova`. Na koncu pa rešujemo sistem $L(Ux^{n+1})=x^{n}$ znotraj metode `inv_lastni`.


## Zgornje Hessenbergove Matrike

Zgornje Hessenbergova matrika je poseben tip kvadratne matrike, pri kateri so vsi elementi pod spodnjo diagonalo enaki nič. Takšna struktura je uporabna, ker izračunavanje LU razcepa zahteva manj korakov. 
Zgornje Hessenbergova matrika se dobi iz Hessenbergovega razcepa določene matrike na različne načine. V tem modulu pa je Hessenbergov razcep implementiran s Hausholderjevimi zrcaljenji.

## Funkcije

### `hessenberg(A::Matrix{T}) -> ZgornjiHessenberg{T}, Matrix{T}`

Izračuna zgornje Hessenbergov razcep matrike `A` in dimenzije `nxn`. Vrne zgornje Hessenbergovo matriko ter matriko prehoda `Q`, ki omogoča transformacijo matrike `A` v zgornje Hessenbergovo matriko.
Householderjeva zrcaljenja potekajo tako:
- iz matrike A vzamemo vektor pod $p$-tim elementom na glavni diagonali (v prvi iteraciji to je prvi element):
```math
\begin{bmatrix}
a & b & \cdots & c & d \\
x & e & \cdots & f & g \\
x & h & \cdots & i & j \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
x & w & \cdots & y & z
\end{bmatrix}
```
- dobimo nek vektor $x$
- ustvarimo vektor $w$, kjer je prva komponenta različnega predznaka kot prva komponenta vektora $x$:
```math
\begin{bmatrix}
\pm \|x\| \\
0 \\
0
\end{bmatrix}
```
- izračunamo vektor $v = w - x$
- ustvarimo matriko $P = \frac{vv^T}{v^Tv}$
- izračunamo matriko $\hat{H} = I - 2P$, ki je ortogonalno simetrična matrika
- vzamemo spodnji desni blok matrike dimenzij $n - p$ in ga vstavimo namesto spodnjega desnega bloka matrike $I$ da dobimo matriko $Q$
- z reševanjem produkta $QAQ$ (ker je matrika ortogonalno simetrična) dobimo ničelne elemente pod spodnjo diagonalo.
- postopek ponovimo za ostale stolpce in na kocu dobimo zgornje Hessenbergovo matriko $H$

#### Argumenti

- `A::Matrix{T}`: Vhodna matrika, ki jo želimo razcepiti.

#### Vrnjeni rezultati

- `ZgornjiHessenberg{T}`: Struktura, ki vsebuje zgornje Hessenbergovo matriko `H`.
- `Matrix{T}`: Matrika prehoda `Q`.

### `lu(A::ZgornjiHessenberg{T}) -> SpodnjaTridiagonalna{T}, Matrix{T}, Matrix{Int}`

Izvede LU razcep nad zgornje Hessenbergovo matriko `A.H`. Vrne razcepljene matrike `L` in `U`, ter matriko permutacij `P`. 
Zaradi ničel pod spodnjo diagonalo je manj operacij. Razcep lahko izračunamo z algoritmom:
```
algorithm MyAlgorithm
    U = copy(H)
    for i = 1: n-1 do:
        l_i = h_i+1_i / u_i_i
        for j = i + 1 : n do
            u_i+1_j = h_i+1_j - l_i*u_i_j
        end for
        append l_i to L
    end for
    return L
end algorithm
```
Zaradi majhnega števila nenačelnih elementov, v vektorju L hranimo samo neničelne elemente. Matrika L je spodnje trikotna in če je obrnljiva je bidiagonalna. Ko na vhodu imamo podano singularno matriko, potem jo z delnim pivotarjenjem preoblikujemo in naredimo zgornje prikazani algoritem. Matriko L hranimo kot vektor neničelnih elementov.

#### Argumenti

- `A::ZgornjiHessenberg{T}`: Struktura s zgornje Hessenbergovo matriko.

#### Vrnjeni rezultati

- `SpodnjaTridiagonalna{T}`: Struktura, ki vsebuje neničelne elemente pod diagonalo matrike `L`.
- `Matrix{T}`: Zgornjetrikotna matrika `U`.
- `Matrix{Int}`: Matrika permutacij `P`.

### `inv_lastni(A::Matrix{T}, l::Real) -> Vector{T}, Real`

Rešuje problem za prirejeno matriko polinoma, ki je podana z matriko `A` in začetnim približkom za lastno vrednost `l`. Vrne lastni vektor in lastno vrednost. 
Metoda implementira algoritem [*Algorithm 18.6: Inverse Iteration to Find Eigenvector of an Upper Hessenberg Matrix*](https://www.sciencedirect.com/topics/mathematics/inverse-power-method):
![](https://github.com/bi4528/ipm/blob/master/IpmZgornjeHessenbergova/ipm_alg.png)

V implementaciji je posebna pozornost namenjena, če je matrika U iz LU razcepa zgornje Hessenbergove matrike singularna. V tem primeru se reševanje sistema $L(Ux^{n+1})=x^{n}$ ne izvaja, ker je `l` že lastna vrednost matrike `A`.

Zaradi implementacije matrike `L` kot vektorja, je definiran poseben operator `\`, ki uporablja dejstvo, da je matrika `L` bidiagonalna. Zaradi tega reševanje `y = L \ xn`, kjer je `xn` lastni vektor poteka kot `y_i+1 = y_i+1 - l_i*x_i`.

#### Argumenti

- `A::Matrix{T}`: Matrika prirejenega polinoma.
- `l::Real`: Začetni približek za lastno vrednost.

#### Vrnjeni rezultati

- `Vector{T}`: Lastni vektor.
- `Real`: Lastna vrednost.

### Primeri delovanja
V skripti `IpmZgornjeHessenbergova/scripts/ipm.jl` si lahko ogledate primere delovanja metod.

## Testiranje in rezultati

Inverzna potenčna metoda `inv_lastni` je testirana na polinomih z realnimi ničlami. Ob testiranju se poljubno ustvari vektor ničel polinoma `roots`, iz katerega se potem ustvari polinom in se koeficienti le-tega shranijo v drug vektor `coeffs`. Tem koeficientom se priredi matrika `B` in na njej se izvede funkcija `inv_lastni` s poljubnim začetnim približkom. Pridobljeni približek lastne vrednosti se potem primerja z elementi vektorja `roots` in če je razlika manjša od $1^{-10}$, se sprejema kot lastna vrednost matrike.
Rezultati kažejo, da za manjše matrike (do recimo `10x10`) uspešno vrača dobre približke lastnih vrednosti, medtem ko za velike matrike z elementi, ki so si zelo blizu lahko vrne slabše približke (recimo `100x100`) in zelo počasi vrača rezultat.

Dokumentacijo je pripravil [Bojan Ilić].
