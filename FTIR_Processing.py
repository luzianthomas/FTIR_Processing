# Processing time domain interferograms recorded on a FT-IR spectrometer

# -----------------
# Imports
# -----------------
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


# -----------------
# Functions
# -----------------
def load_data(filename):
    return np.loadtxt(filename, delimiter=",")
    
def extract_zpd_region(signal, N):
    ix = np.argmax(signal[:,1])
    low = int(ix - N/2)
    high = int(ix + N/2)
    return signal[low:high,1]

def extract_and_average_zpd(signal, N):
    half = len(signal) // 2
    first_half = signal[:half]
    second_half = signal[half:]

    ifg1 = extract_zpd_region(first_half, N)
    ifg2 = extract_zpd_region(second_half, N)

    avg_ifg = (ifg1 + ifg2) / 2
    return avg_ifg, ifg1, ifg2

def compute_spectrum(interferogram, N):
    return np.abs(np.fft.rfft(interferogram, n=N))

def compute_transmittance(sample, background):
    return sample / background

def compute_absorbance(transmittance):
    return -np.log10(transmittance)

def select_region(nu, spectrum, nu_min, nu_max):
    mask = (nu >= nu_min) & (nu <= nu_max)
    return nu[mask], spectrum[mask]

def detect_peaks(absorbance, prominence_factor=0.1):
    threshold = max(absorbance) * prominence_factor
    peaks, _ = find_peaks(absorbance, prominence=threshold)
    return peaks

def print_peak_list(nu, absorbance, peaks_idx, filename=None):
    
    print(f'Peak list:\nWavenumber [cm-1]\tAbsorbance [A]')
    for d in reversed(peaks_idx):
        print(f'{nu[d]:.5f}\t{absorbance[d]:.2f}')
        
    if filename:
        with open(filename, "w") as f:
            f.write("Wavenumber (cm^-1)\tAbsorbance\n")
            for d in reversed(peaks_idx):
                f.write(f'{nu[d]:.5f}\t{absorbance[d]:.2f}\n')

def plot_spectrum(nu, spectrum, peaks=None, ylabel="Intensity", title=None, filename=None):
    plt.figure(figsize=(8,5))
    plt.plot(nu, spectrum, label='Signal')
    if peaks is not None:
        plt.plot(nu[peaks], spectrum[peaks], 'x', color='red', label='Peaks')
    plt.xlabel(r"Wavenumber (cm$^{-1}$)")
    plt.ylabel(ylabel)
    if title:
        plt.title(title)
    plt.gca().invert_xaxis()  # typische Darstellung in FTIR
    plt.tight_layout()
    if filename:
        plt.savefig(filename)
    plt.show()

# -----------------
# Main logic
# -----------------
def main():
    # Parameter
    N = 2**14   # Points for fft
    nu_nyq = 16707.63   # Nyquist frequency of the device (cm-1)
    nu_min = 400   # Minimum frequency (cm-1) for peak picking and plot
    nu_max = 4000   # Maximum frequency (cm-1) for peak picking and plot
    prominence = 0.03  # Prominence factor for prominence based peak detection
    
    # load data from files
    bg = load_data("IRData/BAlk.RIFG.dpt")
    sample = load_data("IRData/BAlk.SIFG.dpt")
    
    # Selection of region around ZPD and averaging
    bg_ifg_avg, _, _ = extract_and_average_zpd(bg, N)
    sample_ifg_avg, _, _ = extract_and_average_zpd(sample, N)

    # fft
    bg_spec = compute_spectrum(bg_ifg_avg, N)
    sample_spec = compute_spectrum(sample_ifg_avg, N)

    # frequency axis
    nu = np.linspace(0, nu_nyq, len(bg_spec))

    # Transmittance and Absorbance
    transmittance = compute_transmittance(sample_spec, bg_spec)
    absorbance = compute_absorbance(transmittance)

    # Selects the relevant region (most commonly 400 to 4000 cm-1)
    nu_region, transmittance_region = select_region(nu, transmittance, nu_min, nu_max)
    nu_region, absorbance_region = select_region(nu, absorbance, nu_min, nu_max)

    # Peak picking
    peaks_idx = detect_peaks(absorbance_region, prominence_factor=prominence)

    # Plots
    plot_spectrum(nu_region, transmittance_region*100, ylabel="Transmittance (%T)", title="Transmission Spectrum")
    plot_spectrum(nu_region, absorbance_region, peaks_idx, ylabel="Absorbance (A)", title="Absorbance Spectrum with Peaks")

    # Print peak list to console (and text file if given)
    print_peak_list(nu_region, absorbance_region, peaks_idx, filename=None)

    
# -----------------
# Entry point
# -----------------
if __name__ == "__main__":
    main()

