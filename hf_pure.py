
"""
Eigenstaendige Hartree-Fock-Implementierung (STO-3G, s- und p-Funktionen)
nur mit numpy und der Python-Standardbibliothek (math.erf fuer die Boys-
Funktion). McMurchie-Davidson-Schema fuer die Gauss-Integrale.
"""
import math
import numpy as np
import itertools
from functools import lru_cache
 
# ------------------------------------------------------------------
# STO-3G Basissatz (exakte Standard-Parameter, Hehre/Stewart/Pople 1969)
# Jede Schale: (drehimpuls 0='s' 1='p', [(exponent, kontraktionskoeff), ...])
# ------------------------------------------------------------------
STO3G = {
    "H": [
        (0, [(3.42525091, 0.15432897), (0.62391373, 0.53532814), (0.1688554, 0.44463454)]),
    ],
    "N": [
        (0, [(99.106169, 0.15432897), (18.052312, 0.53532814), (4.8856602, 0.44463454)]),
        (0, [(3.7804559, -0.09996723), (0.8784966, 0.39951283), (0.2857144, 0.70011547)]),
        (1, [(3.7804559, 0.15591627), (0.8784966, 0.60768372), (0.2857144, 0.39195739)]),
    ],
    "C": [
        (0, [(71.616837, 0.15432897), (13.045096, 0.53532814), (3.5305122, 0.44463454)]),
        (0, [(2.9412494, -0.09996723), (0.6834831, 0.39951283), (0.2222899, 0.70011547)]),
        (1, [(2.9412494, 0.15591627), (0.6834831, 0.60768372), (0.2222899, 0.39195739)]),
    ],
    "O": [
        (0, [(130.70932, 0.15432897), (23.808861, 0.53532814), (6.4436083, 0.44463454)]),
        (0, [(5.0331513, -0.09996723), (1.1695961, 0.39951283), (0.380389, 0.70011547)]),
        (1, [(5.0331513, 0.15591627), (1.1695961, 0.60768372), (0.380389, 0.39195739)]),
    ],
}
ATOM_Z = {"H": 1, "C": 6, "N": 7, "O": 8}
 
CART_P = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]  # px, py, pz
 
 
def primitive_norm(alpha, l, m, n):
    """Normierung einer einzelnen primitiven kartesischen Gauss-Funktion."""
    L = l + m + n
    num = (2 * alpha / math.pi) ** 0.75 * (4 * alpha) ** (L / 2)
    den = math.sqrt(
        math.factorial(2 * l) // (2 ** l * math.factorial(l)) *
        math.factorial(2 * m) // (2 ** m * math.factorial(m)) *
        math.factorial(2 * n) // (2 ** n * math.factorial(n))
        if False else 1
    )
    # (2l-1)!! als math.factorial(2l)/(2^l l!)
    def doublefact_odd(k):
        # (2k-1)!!
        r = 1
        i = 2 * k - 1
        while i > 0:
            r *= i
            i -= 2
        return r
    denom = math.sqrt(doublefact_odd(l) * doublefact_odd(m) * doublefact_odd(n))
    return num / denom
 
 
class Shell:
    def __init__(self, center, ang, prims):
        self.center = np.array(center, dtype=float)
        self.ang = ang  # 0 = s, 1 = p
        # prims: Liste (exponent, kontraktionskoeff) - noch nicht normiert
        self.prims = prims
 
 
def build_shells(atoms):
    """atoms: Liste (Symbol, (x,y,z)) in Angstrom -> Liste Shell-Objekte."""
    ANG2BOHR = 1.8897259886
    shells = []
    for sym, pos in atoms:
        pos_bohr = np.array(pos) * ANG2BOHR
        for ang, prims in STO3G[sym]:
            shells.append(Shell(pos_bohr, ang, prims))
    return shells
 
 
def contracted_norm(ang, prims):
    """Normierungsfaktor je Primitive, dann Gesamtnormierung der Kontraktion,
    damit die Selbstueberlappung der Kontraktion (fuer den ersten kartesischen
    Kanal, z.B. px) exakt 1 ist."""
    l, m, n = (ang, 0, 0)  # reicht fuer s (0,0,0) und einen p-Kanal (1,0,0);
    # Normierungsfaktor ist fuer alle drei p-Kanaele wegen Symmetrie gleich.
    coeffs_norm = []
    for alpha, c in prims:
        coeffs_norm.append(c * primitive_norm(alpha, l, m, n))
    # Selbstueberlappung der kontrahierten (noch unnormierten) Funktion
    s = 0.0
    for (a1, c1), (a2, c2) in itertools.product(zip([p[0] for p in prims], coeffs_norm),
                                                 repeat=2):
        p_ = a1 + a2
        s += c1 * c2 * (math.pi / p_) ** 1.5 * doublefact_ratio(l, m, n, p_)
    norm = 1.0 / math.sqrt(s)
    return [c * norm for c in coeffs_norm]
 
 
def doublefact_ratio(l, m, n, p):
    """Fuer den Selbstueberlappungs-Ausdruck einer kartesischen Gauss-Funktion
    mit gleichem Exponenten auf beiden Seiten (Zentrum=0 relativ)."""
    def doublefact_odd(k):
        r = 1
        i = 2 * k - 1
        while i > 0:
            r *= i
            i -= 2
        return r
    return (doublefact_odd(l) * doublefact_odd(m) * doublefact_odd(n)) / (2 * p) ** (l + m + n)
 
 
# ------------------------------------------------------------------
# McMurchie-Davidson Hermite-Koeffizienten
# ------------------------------------------------------------------
@lru_cache(maxsize=None)
def E(i, j, t, Qx, a, b):
    p = a + b
    q = a * b / p
    if t < 0 or t > i + j:
        return 0.0
    if i == 0 and j == 0 and t == 0:
        return math.exp(-q * Qx * Qx)
    if j == 0:
        return (1 / (2 * p)) * E(i - 1, j, t - 1, Qx, a, b) \
            - (q * Qx / a) * E(i - 1, j, t, Qx, a, b) \
            + (t + 1) * E(i - 1, j, t + 1, Qx, a, b)
    else:
        return (1 / (2 * p)) * E(i, j - 1, t - 1, Qx, a, b) \
            + (q * Qx / b) * E(i, j - 1, t, Qx, a, b) \
            + (t + 1) * E(i, j - 1, t + 1, Qx, a, b)
 
 
def boys(n, x):
    """Boys-Funktion F_n(x), nur mit math (kein scipy)."""
    if x < 1e-12:
        return 1.0 / (2 * n + 1)
    if x < 25:
        # Taylor-Serie fuer hohes n, dann Abwaertsrekursion (stabil)
        nmax = n + 15
        term = 1.0 / (2 * nmax + 1)
        summe = term
        k = 1
        while True:
            term *= x / (nmax + 0.5 + k)
            summe += term
            if term < 1e-17:
                break
            k += 1
        f = math.exp(-x) * summe
        vals = [0.0] * (nmax + 1)
        vals[nmax] = f
        for m in range(nmax, 0, -1):
            vals[m - 1] = (2 * x * vals[m] + math.exp(-x)) / (2 * m - 1)
        return vals[n]
    else:
        # asymptotische Formel ueber die Fehlerfunktion, Aufwaertsrekursion
        f0 = 0.5 * math.sqrt(math.pi / x) * math.erf(math.sqrt(x))
        vals = [0.0] * (n + 1)
        vals[0] = f0
        for m in range(1, n + 1):
            vals[m] = ((2 * m - 1) * vals[m - 1] - math.exp(-x)) / (2 * x)
        return vals[n]
 
 
@lru_cache(maxsize=None)
def hermite_R(t, u, v, n, p, PCx, PCy, PCz):
    T = PCx ** 2 + PCy ** 2 + PCz ** 2
    if t == 0 and u == 0 and v == 0:
        return (-2 * p) ** n * boys(n, p * T)
    if t > 0:
        val = 0.0
        if t > 1:
            val += (t - 1) * hermite_R(t - 2, u, v, n + 1, p, PCx, PCy, PCz)
        val += PCx * hermite_R(t - 1, u, v, n + 1, p, PCx, PCy, PCz)
        return val
    if u > 0:
        val = 0.0
        if u > 1:
            val += (u - 1) * hermite_R(t, u - 2, v, n + 1, p, PCx, PCy, PCz)
        val += PCy * hermite_R(t, u - 1, v, n + 1, p, PCx, PCy, PCz)
        return val
    val = 0.0
    if v > 1:
        val += (v - 1) * hermite_R(t, u, v - 2, n + 1, p, PCx, PCy, PCz)
    val += PCz * hermite_R(t, u, v - 1, n + 1, p, PCx, PCy, PCz)
    return val
 
 
# ------------------------------------------------------------------
# Basisfunktionsliste: (shell_index, kartesischer Kanal (l,m,n))
# ------------------------------------------------------------------
def expand_basis_functions(shells):
    funcs = []  # (shell, (l,m,n), normierte_coeffs)
    for sh in shells:
        if sh.ang == 0:
            norm_c = contracted_norm(0, sh.prims)
            funcs.append((sh, (0, 0, 0), norm_c))
        else:
            norm_c = contracted_norm(1, sh.prims)
            for lmn in CART_P:
                funcs.append((sh, lmn, norm_c))
    return funcs
 
 
def prim_pair_terms(fa, fb):
    """Alle primitiven Paar-Terme (exp_a, ca, exp_b, cb) zweier Basisfunktionen."""
    sh_a, lmn_a, ca_list = fa
    sh_b, lmn_b, cb_list = fb
    out = []
    for (alpha, _), ca in zip(sh_a.prims, ca_list):
        for (beta, _), cb in zip(sh_b.prims, cb_list):
            out.append((alpha, ca, beta, cb))
    return out
 
 
def overlap_prim(alpha, A, la, beta, B, lb):
    p = alpha + beta
    AB = A - B
    Sx = E(la[0], lb[0], 0, AB[0], alpha, beta)
    Sy = E(la[1], lb[1], 0, AB[1], alpha, beta)
    Sz = E(la[2], lb[2], 0, AB[2], alpha, beta)
    return Sx * Sy * Sz * (math.pi / p) ** 1.5
 
 
def kinetic_prim(alpha, A, la, beta, B, lb):
    p = alpha + beta
    AB = A - B
 
    def S1(i, j, dim):
        return E(i, j, 0, AB[dim], alpha, beta)
 
    def T1(i, j, dim):
        term = beta * (2 * j + 1) * S1(i, j, dim)
        term -= 2 * beta ** 2 * S1(i, j + 2, dim)
        if j >= 2:
            term -= 0.5 * j * (j - 1) * S1(i, j - 2, dim)
        return term
 
    pre = (math.pi / p) ** 1.5
    Sx, Sy, Sz = S1(la[0], lb[0], 0), S1(la[1], lb[1], 1), S1(la[2], lb[2], 2)
    Tx, Ty, Tz = T1(la[0], lb[0], 0), T1(la[1], lb[1], 1), T1(la[2], lb[2], 2)
    return pre * (Tx * Sy * Sz + Sx * Ty * Sz + Sx * Sy * Tz)
 
 
def nuclear_prim(alpha, A, la, beta, B, lb, C):
    p = alpha + beta
    P = (alpha * A + beta * B) / p
    AB = A - B
    PC = P - C
    val = 0.0
    for t in range(la[0] + lb[0] + 1):
        Ex = E(la[0], lb[0], t, AB[0], alpha, beta)
        for u in range(la[1] + lb[1] + 1):
            Ey = E(la[1], lb[1], u, AB[1], alpha, beta)
            for v in range(la[2] + lb[2] + 1):
                Ez = E(la[2], lb[2], v, AB[2], alpha, beta)
                val += Ex * Ey * Ez * hermite_R(t, u, v, 0, p, PC[0], PC[1], PC[2])
    return val * 2 * math.pi / p
 
 
def eri_prim(alpha, A, la, beta, B, lb, gamma, C, lc, delta, D, ld):
    p = alpha + beta
    q = gamma + delta
    P = (alpha * A + beta * B) / p
    Q = (gamma * C + delta * D) / q
    AB = A - B
    CD = C - D
    PQ = P - Q
    alpha_tot = p * q / (p + q)
    val = 0.0
    for t1 in range(la[0] + lb[0] + 1):
        Ex1 = E(la[0], lb[0], t1, AB[0], alpha, beta)
        for u1 in range(la[1] + lb[1] + 1):
            Ey1 = E(la[1], lb[1], u1, AB[1], alpha, beta)
            for v1 in range(la[2] + lb[2] + 1):
                Ez1 = E(la[2], lb[2], v1, AB[2], alpha, beta)
                E1 = Ex1 * Ey1 * Ez1
                if E1 == 0.0:
                    continue
                for t2 in range(lc[0] + ld[0] + 1):
                    Ex2 = E(lc[0], ld[0], t2, CD[0], gamma, delta)
                    for u2 in range(lc[1] + ld[1] + 1):
                        Ey2 = E(lc[1], ld[1], u2, CD[1], gamma, delta)
                        for v2 in range(lc[2] + ld[2] + 1):
                            Ez2 = E(lc[2], ld[2], v2, CD[2], gamma, delta)
                            E2 = Ex2 * Ey2 * Ez2
                            if E2 == 0.0:
                                continue
                            sign = (-1) ** (t2 + u2 + v2)
                            R = hermite_R(t1 + t2, u1 + u2, v1 + v2, 0,
                                          alpha_tot, PQ[0], PQ[1], PQ[2])
                            val += E1 * E2 * sign * R
    return val * 2 * math.pi ** 2.5 / (p * q * math.sqrt(p + q))
 
 
def build_integrals(atoms):
    shells = build_shells(atoms)
    funcs = expand_basis_functions(shells)
    n = len(funcs)
    S = np.zeros((n, n))
    T = np.zeros((n, n))
    V = np.zeros((n, n))
    charges = [(np.array(pos) * 1.8897259886, ATOM_Z[sym]) for sym, pos in atoms]
 
    for i, fa in enumerate(funcs):
        for j, fb in enumerate(funcs):
            if j < i:
                continue
            sh_a, la, ca = fa
            sh_b, lb, cb = fb
            s_val = t_val = v_val = 0.0
            for alpha, ca_i, beta, cb_j in prim_pair_terms(fa, fb):
                cc = ca_i * cb_j
                s_val += cc * overlap_prim(alpha, sh_a.center, la, beta, sh_b.center, lb)
                t_val += cc * kinetic_prim(alpha, sh_a.center, la, beta, sh_b.center, lb)
                for C, Z in charges:
                    v_val -= cc * Z * nuclear_prim(alpha, sh_a.center, la, beta, sh_b.center, lb, C)
            S[i, j] = S[j, i] = s_val
            T[i, j] = T[j, i] = t_val
            V[i, j] = V[j, i] = v_val
 
    ERI = np.zeros((n, n, n, n))
    for i, fa in enumerate(funcs):
        for j, fb in enumerate(funcs):
            if j > i:
                continue
            for k, fc in enumerate(funcs):
                for l, fd in enumerate(funcs):
                    if l > k or (k * n + l) > (i * n + j):
                        continue
                    sh_a, la, ca = fa
                    sh_b, lb, cb = fb
                    sh_c, lc_, cc_ = fc
                    sh_d, ld, cd = fd
                    val = 0.0
                    for a1, ca1, b1, cb1 in prim_pair_terms(fa, fb):
                        for a2, cc2, b2, cd2 in prim_pair_terms(fc, fd):
                            coeff = ca1 * cb1 * cc2 * cd2
                            val += coeff * eri_prim(a1, sh_a.center, la, b1, sh_b.center, lb,
                                                     a2, sh_c.center, lc_, b2, sh_d.center, ld)
                    for (p1, q1, r1, s1) in {(i, j, k, l), (j, i, k, l), (i, j, l, k), (j, i, l, k),
                                              (k, l, i, j), (l, k, i, j), (k, l, j, i), (l, k, j, i)}:
                        ERI[p1, q1, r1, s1] = val
    return S, T, V, ERI, n
 
 
def _kernabstossung(atoms):
    E_nuc = 0.0
    coords = [(np.array(pos) * 1.8897259886, ATOM_Z[sym]) for sym, pos in atoms]
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            Ci, Zi = coords[i]
            Cj, Zj = coords[j]
            E_nuc += Zi * Zj / np.linalg.norm(Ci - Cj)
    return E_nuc
 
 
def _scf_lauf(S, Hcore, ERI, nbf, n_alpha, n_beta, mix, max_iter, tol):
    """Ein SCF-Durchlauf mit fester Daempfung. Gibt (E_elec, konvergiert)."""
    eigval, eigvec = np.linalg.eigh(S)
    X = eigvec @ np.diag(eigval ** -0.5) @ eigvec.T
 
    def build_density(C, nocc):
        Cocc = C[:, :nocc]
        return Cocc @ Cocc.T
 
    Pa = np.zeros((nbf, nbf))
    Pb = np.zeros((nbf, nbf))
    E_old = 0.0
    E_elec = 0.0
    letzte = []
 
    for it in range(max_iter):
        J = np.einsum('pqrs,rs->pq', ERI, Pa + Pb)
        Ka = np.einsum('prqs,rs->pq', ERI, Pa)
        Kb = np.einsum('prqs,rs->pq', ERI, Pb)
        Fa = Hcore + J - Ka
        Fb = Hcore + J - Kb
 
        ea, Ca_p = np.linalg.eigh(X.T @ Fa @ X)
        eb, Cb_p = np.linalg.eigh(X.T @ Fb @ X)
        Ca, Cb = X @ Ca_p, X @ Cb_p
 
        Pa_ziel = build_density(Ca, n_alpha)
        Pb_ziel = build_density(Cb, n_beta)
 
        E_elec = 0.5 * np.sum((Pa_ziel + Pb_ziel) * Hcore) \
            + 0.5 * np.sum(Pa_ziel * Fa) + 0.5 * np.sum(Pb_ziel * Fb)
        dP = np.linalg.norm(Pa_ziel - Pa) + np.linalg.norm(Pb_ziel - Pb)
 
        Pa = mix * Pa_ziel + (1 - mix) * Pa
        Pb = mix * Pb_ziel + (1 - mix) * Pb
 
        letzte.append(E_elec)
        if abs(E_elec - E_old) < tol and dP < 1e-7 and it > 3:
            return E_elec, True
        E_old = E_elec
 
    # Nicht innerhalb max_iter konvergiert - als letzte Zahlen unruhig?
    stabil = len(letzte) >= 8 and np.std(letzte[-8:]) < 1e-6
    return E_elec, stabil
 
 
def scf(atoms, n_electrons, spin=0, max_iter=150, tol=1e-9):
    S, T, V, ERI, nbf = build_integrals(atoms)
    Hcore = T + V
    n_alpha = (n_electrons + spin) // 2
    n_beta = n_electrons - n_alpha
 
    # Erst vorsichtig (starke Daempfung), dann ggf. mit mehr Schritten -
    # das verhindert, dass die Loesung bei schwierigen Geometrien in ein
    # falsches (oszillierendes) Fixpunkt-Paar einrastet.
    for mix, iters in [(0.3, max_iter), (0.15, max_iter * 3)]:
        E_elec, ok = _scf_lauf(S, Hcore, ERI, nbf, n_alpha, n_beta,
                                mix, iters, tol)
        if ok:
            break
 
    return E_elec + _kernabstossung(atoms)
 
 
def energie(atoms, spin=0):
    """Bequemer Einstiegspunkt: Gesamtenergie eines Molekuels (STO-3G HF).
 
    atoms: Liste (Elementsymbol, (x,y,z)) in Angstrom, z.B.
           [("H", (0,0,0)), ("H", (0,0,0.74))]
    spin:  Anzahl ungepaarter Elektronen (0 fuer die meisten Molekuele,
           1 fuer ein einzelnes H-Atom, 2 fuer C- oder O-Atome).
    Gibt die Energie in Hartree zurueck.
    """
    n_elec = sum(ATOM_Z[sym] for sym, _ in atoms)
    return scf(atoms, n_electrons=n_elec, spin=spin)
