import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt


import pyswarms as ps
from pyswarms.single.global_best import GlobalBestPSO
import pycbc
import pycbc.waveform
import pycbc.psd


import GWFish.modules as gw
import GWFish.modules.constants as cst


def get_duration_between_two_freqs_seconds(m1, m2, freq1, freq2):
    """
    Calculate the time duration between two gravitational wave frequencies
    for a binary system with given component masses.

    Parameters:
    m1 (float): Mass of the first object in solar masses.
    m2 (float): Mass of the second object in solar masses.
    freq1 (float): The initial frequency in Hz.
    freq2 (float): The final frequency in Hz.

    Returns:
    float: The absolute time duration between the two frequencies in seconds.
    """
    Msol = 1.9885e30
    c = 299792458.
    G = 6.674e-11

    m1 *= Msol
    m2 *= Msol

    M = m1 + m2
    mu = m1 * m2 / M

    Mc = G * mu ** 0.6 * M ** 0.4 / c ** 3

    time = -5. / (256. * np.pi ** (8 / 3)) / Mc ** (5 / 3) / freq1 ** (8 / 3)
    time_2 = -5. / (256. * np.pi ** (8 / 3)) / Mc ** (5 / 3) / freq2 ** (8 / 3)

    return np.abs(time - time_2)

def fd_gwfish_output_format(hfp, hfc):
    hfp = hfp[:, np.newaxis]
    hfc = hfc[:, np.newaxis]
    polarizations = np.hstack((hfp, hfc))
    return polarizations

def t_of_f_PN(parameters, frequencyvector):
    local_params = parameters.copy()
    M1 = local_params['mass_1'] * cst.Msol
    M2 = local_params['mass_2'] * cst.Msol
    M = M1 + M2
    mu = M1 * M2 / M
    Mc = cst.G * mu ** 0.6 * M ** 0.4 / cst.c ** 3
    t_of_f = -5. / (256. * np.pi ** (8 / 3)) / Mc ** (5 / 3) / frequencyvector ** (8 / 3)
    return t_of_f + local_params['geocent_time']

def pycbc_match(injection, template, **kwargs):
    f_low = kwargs['f_low']
    f_high = kwargs['f_high']
    psd = pycbc.psd.from_txt(kwargs['psd_path'], len(injection), kwargs['delta_f'], f_low, is_asd_file=False) / kwargs['number_of_stations']
    return pycbc.filter.matchedfilter.match(injection, template, psd=psd, low_frequency_cutoff=f_low, high_frequency_cutoff=f_high)[0]

def match(channel_1, channel_2, template, **kwargs):
    match_1 = pycbc_match(channel_1, template, **kwargs)
    match_2 = pycbc_match(channel_2, template, **kwargs)
    return np.sqrt(match_1 ** 2 + match_2 ** 2) / np.sqrt(2)

def gen_template(max_prms, **temp_kwargs):
    m1, m2 = max_prms
    hp, _ = pycbc.waveform.get_fd_waveform(
        approximant=temp_kwargs['approximant'],
        mass1=m1,
        mass2=m2,
        spin1x=0,
        spin1y=0,
        spin1z=0,
        spin2x=0,
        spin2y=0,
        spin2z=0,
        inclination=temp_kwargs['theta_jn'],
        coa_phase=temp_kwargs['phase'],
        distance=temp_kwargs['luminosity_distance'],
        eccentricity=temp_kwargs['eccentricity'],
        delta_f=freq_prms['delta_f'],
        f_lower=freq_prms['f_low_injection'],
        f_final=freq_prms['f_high'],
        f_ref=freq_prms['f_ref']
    )
    return hp

def is_valid_mass(m1, m2):
    return m1 >= 3.5 and m2 >= 3.5

def objective_match(x, injection_1, injection_2, **kwargs):
    m1, m2 = gw.auxiliary.from_mChirp_eta_to_m1_m2(x[0], x[1])
    
    if not is_valid_mass(m1, m2):
        return 1e4

    template = gen_template((m1, m2), **kwargs)
    return -1 * match(injection_1, injection_2, template, **kwargs)

def objective_wrapper(x, injection_1, injection_2, **kwargs):
    return np.array([objective_match(particle, injection_1, injection_2, **kwargs) for particle in x])

# Parameter Dictionaries
injection_prms = {
    'chirp_mass': 400,
    'eta': 0.249,
    'a_1': 0,
    'a_2': 0,
    'tilt_1': 0,
    'tilt_2': 0,
    'phi_12': 0,
    'phi_jl': 0,
    'theta_jn': 2.331013,
    'luminosity_distance': 402.237232,
    'ra': 1.487217,
    'dec': -1.25711,
    'psi': 1.7952,
    'phase': 1.469814,
    'geocent_time': 1126259462.4116447,
    'redshift': 0,
    'eccentricity': 0.1
}

injection_prms['mass_1'], injection_prms['mass_2'] = gw.auxiliary.from_mChirp_eta_to_m1_m2(injection_prms['chirp_mass'], injection_prms['eta'])

freq_prms = {
    'f_low_injection': 0.08,
    'f_low': 0.1,
    'f_high': 3,
    'f_ref': 0.1
}

freq_prms['delta_f'] = 1 / get_duration_between_two_freqs_seconds(injection_prms['mass_1'], injection_prms['mass_2'], freq_prms['f_low'], freq_prms['f_high'])

inj_waveform_prms = {
    'approximant': 'TaylorF2Ecc'
}

match_kwargs = {
    **freq_prms,
    'psd_path': '/home/anemmani/softwares/anaconda3/envs/gwfish/lib/python3.10/site-packages/GWFish/detector_psd/LGWA_Si_psd.txt',
    'number_of_stations': 4
}

template_prms = {
    'approximant': 'TaylorF2Ecc',
    'theta_jn': 2.331013,
    'luminosity_distance': 402.237232,
    'ra': 1.487217,
    'dec': -1.25711,
    'psi': 1.7952,
    'phase': 1.469814,
    'geocent_time': 1126259462.4116447,
    'redshift': 0,
    'eccentricity': 0
}

objective_kwargs = {**match_kwargs, **template_prms}

# Injected Signal
hp, hc = pycbc.waveform.get_fd_waveform(
    approximant=inj_waveform_prms['approximant'],
    mass1=injection_prms['mass_1'],
    mass2=injection_prms['mass_2'],
    spin1x=0,
    spin1y=0,
    spin1z=0,
    spin2x=0,
    spin2y=0,
    spin2z=0,
    inclination=injection_prms['theta_jn'],
    coa_phase=injection_prms['phase'],
    distance=injection_prms['luminosity_distance'],
    eccentricity=injection_prms['eccentricity'],
    delta_f=freq_prms['delta_f'],
    f_lower=freq_prms['f_low_injection'],
    f_final=freq_prms['f_high'],
    f_ref=freq_prms['f_ref']
)

polarizations = fd_gwfish_output_format(np.array(hp), np.array(hc))

detector = gw.detection.Detector('LGWA', override_frequencyvector=np.array(hp.sample_frequencies))

timevector = t_of_f_PN(injection_prms, np.array(hp.sample_frequencies))

injection_signal = gw.detection.projection_moon(injection_prms, detector, polarizations, timevector)

injection_channel_1 = pycbc.types.frequencyseries.FrequencySeries(injection_signal[:, 0], delta_f=hp.delta_f)

injection_channel_2 = pycbc.types.frequencyseries.FrequencySeries(injection_signal[:, 1], delta_f=hp.delta_f)

# PSO
x_min = np.array([100, 0.1])
x_max = np.array([700, 0.25])

sigma_mchirp = 10
sigma_eta = 0.1

n_particles = 100

chirp_mass_guess = stats.truncnorm(
    (x_min[0] - injection_prms['chirp_mass']) / sigma_mchirp,
    (x_max[0] - injection_prms['chirp_mass']) / sigma_mchirp,
    loc=injection_prms['chirp_mass'], scale=sigma_mchirp).rvs(n_particles)

eta_guess = stats.truncnorm(
    (x_min[1] - injection_prms['eta']) / sigma_eta,
    (x_max[1] - injection_prms['eta']) / sigma_eta,
    loc=injection_prms['eta'], scale=sigma_eta).rvs(n_particles)

init_pos = np.array([chirp_mass_guess, eta_guess]).T

PSO_OPTIONS = {'c1': 0.3, 'c2': 0.9, 'w': 0.9}

bounds = (x_min, x_max)
optimizer = GlobalBestPSO(n_particles=n_particles, dimensions=2, options=PSO_OPTIONS, bounds=bounds, init_pos=init_pos)

negative_match, rec_prms = optimizer.optimize(objective_wrapper, iters=20, injection_1 = injection_channel_1, injection_2 = injection_channel_2, **objective_kwargs)

print('The Fitting factor match is: ', -1*negative_match)
print('The recovered parameters are: ', rec_prms)

