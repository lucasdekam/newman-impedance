import numpy as np
from scipy.special import factorial
from scipy.special import legendre
from scipy.integrate import solve_bvp

N_INTEGRATION_PTS = 1000  # Number of points for numerical integration of Lagrange function products. 1000 seemed fine from plot observation


def integrate_single_legendre_polynomial(n: int) -> float:
    """
    Solve the integral between 0 and 1 of eta * P_2n (eta)
    """
    eta = np.linspace(0, 1, N_INTEGRATION_PTS)
    y = eta * legendre(2 * n)(eta)
    return np.trapz(y, eta)


def integrate_double_legendre_polynomial(n: int, m: int) -> float:
    """
    Solve the integral between 0 and 1 of eta * P_2n (eta) * P_2m (eta)
    """
    eta = np.linspace(0, 1, N_INTEGRATION_PTS)
    y = eta * legendre(2 * n)(eta) * legendre(2 * m)(eta)
    return np.trapz(y, eta)


def get_vector(nmax: int):
    """
    Get the vector of integrals from 0 to 1 of eta * P_2n (eta) for
    n between 0 and nmax
    """
    vector = np.zeros(nmax)
    for i in range(nmax):
        vector[i] = integrate_single_legendre_polynomial(i)
    return vector


def get_product_matrix(nmax: int):
    """
    Get the matrix of integrals from 0 to 1  eta * P_2n (eta) * P_2m (eta)
    for n and m between 0 and nmax
    """
    matrix = np.zeros((nmax, nmax))
    for i in range(nmax):
        for j in range(nmax):
            if j >= i:
                matrix[i, j] = integrate_double_legendre_polynomial(i, j)
                matrix[j, i] = matrix[i, j]

    return matrix


def get_prefactor_matrix(nmax: int, omega: float, faradaic: float = 0):
    """
    Get the diagonal matrix with entries 1/(j Omega C_m) where C_m is the prefactor
    C_m = -(4m+1)/M'_2m(0)
    where M'_2m(0) is the derivative of a Legendre function evaluated at zero. An expression
    can be found in one of Newman's papers: doi.org/10.1149/1.2423795
    """

    def get_prefactor(n) -> complex:
        m_prime = -2 / np.pi * (2**n * factorial(n)) ** 4 / factorial(2 * n) ** 2
        return -(4 * n + 1) / m_prime

    diag = np.array(
        [1 / ((1j * omega + faradaic) * get_prefactor(i)) for i in range(nmax)]
    )

    return np.diag(diag)


def calculate_b0(nmax: int, omega: np.ndarray, faradaic: float = 0):
    """
    Calculate the coefficients B by solving the matrix-vector problem. We only need the first one,
    so return B0 for each value of omega.
    """

    def solve(nmax: int, omega_value):
        return np.linalg.solve(
            get_prefactor_matrix(nmax, omega_value, faradaic)
            + get_product_matrix(nmax),
            get_vector(nmax),
        )[0]

    return np.array([solve(nmax, x) for x in omega])


def calculate_impedance(
    nmax: int,
    freq_list: np.ndarray,
    capacitance: float,
    conductivity: float,
    r0: float,
    faradaic: float = 0,
):
    """
    Calculate the impedance Z ~ 1/B0, with realistic frequencies and other parameters.
    """
    omega_list = 2 * np.pi * freq_list * capacitance * r0 / conductivity
    b0 = calculate_b0(nmax, omega_list, faradaic)
    return 1 / (4 * r0 * conductivity * b0)


def solve_leg(n, x_axis: np.ndarray) -> np.ndarray:
    """
    Solve the Legendre-like differential equation
    d/dx [(1+x^2) dy/dx] = n(n+1) y
    which yields the Legendre-like polynomials M_n
    """

    def rhs(x, y):
        ret = np.zeros([2, len(x)])
        ret[0, :] = 1 / (1 + x**2) * y[1, :]
        ret[1, :] = n * (n + 1) * y[0, :]
        return ret

    def bc(ya, yb):
        return np.array([ya[0] - 1, yb[0] - 0])

    sol = solve_bvp(rhs, bc, x=x_axis, y=np.random.randn(2, len(x_axis)), tol=1e-6)
    return sol.y[0, :]


def calculate_u(nmax: int, omega: float, eta: np.ndarray, xi: np.ndarray):
    """
    Calculate the complex potential field U

    Input:
    nmax, omega parameters as usual
    eta, xi: 1D arrays for the coordinates eta and xi
    """
    b_vec = np.linalg.solve(
        get_prefactor_matrix(nmax, omega) + get_product_matrix(nmax), get_vector(nmax)
    )

    u_grid = np.zeros([len(xi), len(eta)], dtype="complex_")
    for n in range(nmax):
        p2n = legendre(2 * n)(eta)
        m2n = solve_leg(2 * n, xi)
        x, y = np.meshgrid(p2n, m2n)
        u_grid += x * y * b_vec[n]

    return u_grid
