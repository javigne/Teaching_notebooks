import numpy as np
import sympy
import matplotlib.pyplot as plt

# DO NOT CHANGE

def get_compton_crosssection(energy_keV, atomic_number, density):
    """reads the nist file, interpolates to energy = energy_keV and return the associated compton cross section
    energy_keV = energy [keV]
    atomic number = Z value
    density = density in [g/cm³]
    returns: cross section [1/cm]
    """
    if energy_keV< 1:
        print('Energy should be larger or equal to 1 keV')
        return 0
    # read in NIST DATA
    input = np.loadtxt('CrossSection_Compton_Z%i.txt' %atomic_number, skiprows =3)
    # set energy in MeV
    energy_MeV = energy_keV/1000
    # interpolate
    cross_section = np.interp(energy_MeV, input[:,0], input[:,1])*density
    return cross_section

def distance_travelled_before_interaction(mu):
    """
    input: cross section for compton scatter [1/cm]
    returns: sampled distance travelled before interaction [cm]
    """
    xi = np.random.random(1)  # one random value between [0,1]
    s = -1/mu * np.log(1-xi)
    return s


# DO NOT CHANGE
''' Kappa
'''


def get_pdf(energy_keV):
    """ the pdf in a function
    """
    mec2 = 5.11e5 # rest energy electron eV
    kappa = lambda E: E * 1000/mec2  # a function that calculates kappa from an incoming energy E in keV
    tau_s = sympy.symbols('tau', real=True, positive=True)  # ratio of the scattered and incoming energy Ec/E
    k_s = sympy.symbols('kappa', real=True, positive=True)  # ratio of the incoming energy and the electron rest energ
    pdf_symp = 1 / tau_s ** 2 + (k_s ** 2 - 2 * k_s - 2) / tau_s + (2 * k_s + 1) + k_s ** 2 * tau_s
    pdf_symp_kfix = pdf_symp.subs(k_s, kappa(energy_keV))  # subsitute 'k_s' by a value
    min_tau = 1 / (1 + 2 * kappa(energy_keV))
    max_tau = 1
    max_value = pdf_symp_kfix.subs(tau_s, min_tau)
    pdf_symp_kfix = pdf_symp_kfix / max_value
    pdf = sympy.lambdify([tau_s], pdf_symp_kfix)
    return pdf, min_tau, max_tau


def rejection_sampling(min_x, max_x, fl_x, print_val=False):
    """
    input
    return: [boolean, sample]
    """
    [xi1, xi2] = np.random.random(2)
    x = min_x + (max_x - min_x) * xi1
    if print_val:
        print('x=', x, " and xi2 =", xi2)
    if xi2 <= fl_x(x):
        return True, x
    else:
        return False, x


def get_phi():
    xi = np.random.random(1)
    phi = np.pi * 2 * xi
    return phi


def get_theta(tau, k):
    cos_theta = 1 - (1 - tau) / (k * tau)
    return np.arccos(cos_theta)


def get_direction(phi, theta):
    de = np.zeros(3)
    de[0] = np.sin(theta) * np.cos(phi)
    de[1] = np.sin(theta) * np.sin(phi)
    de[2] = np.cos(theta)
    return de


def get_detector_position_co(co_i, d, z_detector=10, xs_detector=100, ys_detector=100, print_val=False):
    distance_traveled = (z_detector - co_i[2]) / d[2]
    if not distance_traveled > 0:
        return None
    else:
        co = distance_traveled * d + co_i
        if abs(co[0]) < xs_detector / 2 and abs(co[1]) < ys_detector / 2:
            return co
        else:
            if print_val:
                print('out of bounce detector: (', co[0], ',', co[1], ')')
            return None


class Image:
    def __init__(self, xs_detector, ys_detector, pixelsize, print_val=False):
        self.pixelsize = pixelsize
        self.n_pix_x = int(xs_detector // pixelsize)
        self.n_pix_y = int(ys_detector // pixelsize)
        self.image = np.zeros([self.n_pix_x, self.n_pix_y])
        self.xs_detector = pixelsize * self.n_pix_x
        self.ys_detector = pixelsize * self.n_pix_y
        self.shift_co = np.array([self.xs_detector / 2, self.ys_detector / 2, 0])
        self.print = print_val

    def add_count(self, co):
        co = co + self.shift_co
        pix_idx = int(co[0] // self.pixelsize)
        pix_idy = int(co[1] // self.pixelsize)
        self.image[pix_idx, pix_idy] += 1
        if self.print:
            print('co: ', co[0], ', ', co[1], "added at [", pix_idx, ',', pix_idy, ']')
        return
