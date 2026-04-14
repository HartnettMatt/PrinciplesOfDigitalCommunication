# %% [markdown]
# # Project 1 - Pulse Design (Problem 2.4.A.1)
# BPSK pulse that fits GSM mask, ISI < 0.25dB, max symbol rate

# %%
# Constants and Setup

import numpy as np
from scipy.signal.windows import tukey
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import os

plt.matplotlib.rc("grid", linestyle="dotted")

_script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(_script_dir, "output")
os.makedirs(output_dir, exist_ok=True)

t_u = 40e-6  # max pulse duration
isi_limit_dB = 0.25

# GSM mask from the MATLAB file
freqm_kHz = np.array(
    [-800, -600, -400, -250, -200, -100, 0.0, 100, 200, 250, 400, 600, 800]
)
mask_dB = np.array([-70, -70, -60, -33, -30, 0, 0.0, 0, -30, -33, -60, -70, -70])
freqm_Hz = freqm_kHz * 1e3
mask_interp = interp1d(
    freqm_Hz, mask_dB, kind="linear", bounds_error=False, fill_value=-70.0
)

fs = 10e6  # 10 MHz sample rate
dt = 1.0 / fs
N_pulse = int(t_u * fs)  # 400 samples

t = np.linspace(-t_u / 2, t_u / 2, N_pulse, endpoint=False)

N_fft = 2**14
f = np.fft.fftshift(np.fft.fftfreq(N_fft, d=dt))

# %%
# RRC Pulse (Prob 2.4.B.1)


def rrc_pulse(t_norm, alpha):
    """RRC pulse, normalized domain (t_b=1). Handles singularities at t=0 and t=1/(4a)."""
    u = np.zeros_like(t_norm, dtype=float)

    # t=0 case
    m0 = np.abs(t_norm) < 1e-12
    u[m0] = (1 - alpha) + 4 * alpha / np.pi

    # t = +/-1/(4*alpha) singularity
    if alpha > 1e-12:
        sing = 1.0 / (4 * alpha)
        ms = np.abs(np.abs(t_norm) - sing) < 1e-8
        u[ms] = (alpha / np.sqrt(2)) * (
            (1 - 2 / np.pi) * np.cos(np.pi / (4 * alpha))
            + (1 + 2 / np.pi) * np.sin(np.pi / (4 * alpha))
        )
    else:
        ms = np.zeros_like(t_norm, dtype=bool)

    # everything else
    mg = ~m0 & ~ms
    tg = t_norm[mg]
    num = (1 - alpha) * np.sinc(tg * (1 - alpha)) + (4 * alpha / np.pi) * np.cos(
        np.pi * tg * (1 + alpha)
    )
    den = 1 - (4 * alpha * tg) ** 2
    u[mg] = num / den

    return u


# %%
# Windowed pulse construction (Prob 2.4.3.4)


def build_pulse(t, t_b, alpha, beta):
    """RRC pulse windowed with Tukey, scaled by 1/sqrt(t_b)."""
    t_norm = t / t_b
    u0 = rrc_pulse(t_norm, alpha)
    u0 = u0 / np.sqrt(t_b)
    w = tukey(len(t), alpha=beta)  # beta = taper fraction
    return u0 * w


# %%
# Constraint checks


def check_mask(u):
    """Check if pulse spectrum fits under GSM mask. Returns (pass, spectrum_dB)."""
    U = np.fft.fftshift(np.fft.fft(u, n=N_fft)) * dt
    Upow = np.abs(U) ** 2
    Upow_max = np.max(Upow)
    spec_dB = 10 * np.log10(Upow / Upow_max + 1e-30)
    mask_vals = mask_interp(f)
    return np.all(spec_dB <= mask_vals), spec_dB


def autocorr_lag(u, lag_samp):
    """r_u(tau) at a given lag in samples."""
    N = len(u)
    lag = abs(lag_samp)
    if lag >= N:
        return 0.0
    if lag_samp >= 0:
        return np.sum(u[lag:] * np.conj(u[: N - lag])) * dt
    else:
        return np.sum(u[: N - lag] * np.conj(u[lag:])) * dt


def check_isi(u, t_b):
    """Worst case ISI degradation (Prob 2.4.3.5). Returns (pass, deg_dB)."""
    eu = np.sum(np.abs(u) ** 2) * dt
    samp_per_tb = t_b / dt
    Lmax = int(np.ceil(len(u) * dt / t_b))

    isi_sum = 0.0
    for l1 in range(1, Lmax + 1):
        lag = int(round(l1 * samp_per_tb))
        # both positive and negative offsets
        rp = autocorr_lag(u, lag)
        rn = autocorr_lag(u, -lag)
        isi_sum += abs(np.real(rp)) / eu + abs(np.real(rn)) / eu

    if isi_sum >= 1.0:
        return False, np.inf
    deg = -20 * np.log10(1 - isi_sum)
    return deg <= isi_limit_dB, deg


# %%
# Grid search over (alpha, beta, t_b)

alphas = np.linspace(0.05, 1.0, 30)
betas = np.linspace(0.01, 1.0, 30)
tbs = np.linspace(2e-6, 40e-6, 150)

best_tb = np.inf
best_alpha = None
best_beta = None
best_pulse = None
best_deg = None

print("Coarse Searching...")
total = len(alphas) * len(betas)
done = 0

for a in alphas:
    for b in betas:
        for tb in tbs:
            u = build_pulse(t, tb, a, b)
            ok, _ = check_mask(u)
            if not ok:
                continue
            ok, deg = check_isi(u, tb)
            if not ok:
                continue
            if tb < best_tb:
                best_tb = tb
                best_alpha = a
                best_beta = b
                best_pulse = u.copy()
                best_deg = deg
            break  # smallest passing tb for this (a,b)

print("\n---- Coarse Search Results: ----")
print(f"alpha = {best_alpha:.4f}")
print(f"beta  = {best_beta:.4f}")
print(f"t_b   = {best_tb * 1e6:.2f} us")
print(f"Max symbol rate = {1 / best_tb:.0f} symbols/s")
print(f"ISI degradation = {best_deg:.4f} dB")
# %%
# Fine search around best coarse result

print("Fine searching...")
da = (alphas[1] - alphas[0])
af = np.linspace((best_alpha-da), (best_alpha+da), 60)
db = (betas[1] - betas[0])
bf = np.linspace((best_beta - db), (best_beta + db), 60)
dtb = (tbs[1] - tbs[0]) 
tbf = np.linspace((best_tb - dtb), (best_tb + dtb), 300)

for a in af:
    for b in bf:
        for tb in tbf:
            u = build_pulse(t, tb, a, b)
            ok, _ = check_mask(u)
            if not ok:
                continue
            ok, deg = check_isi(u, tb)
            if not ok:
                continue
            if tb < best_tb:
                best_tb = tb
                best_alpha = a
                best_beta = b
                best_pulse = u.copy()
                best_deg = deg
            break

print("\n---- Fine Search Results: ----")
print(f"alpha = {best_alpha:.4f}")
print(f"beta  = {best_beta:.4f}")
print(f"t_b   = {best_tb * 1e6:.2f} us")
print(f"Max symbol rate = {1 / best_tb:.0f} symbols/s")
print(f"ISI degradation = {best_deg:.4f} dB")


# %%
# Plot: u(t) time domain

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(t * 1e6, best_pulse, "-b", linewidth=1.5)
ax.set_xlabel("Time (us)")
ax.set_ylabel("Amplitude")
ax.set_title(
    f"Optimal Pulse u(t) - a={best_alpha:.3f}, b={best_beta:.3f}, tb={best_tb * 1e6:.2f}us"
)
ax.grid(True)
fig.savefig(os.path.join(output_dir, "pulse_time_domain.png"), dpi=150)
plt.show()
print(f"\nMaximum achievable symbol rate: {1 / best_tb:.0f} symbols/s")

# %%
# Plot: u(t) and r_u(tau)

tau_range = np.linspace(-t_u, t_u, 2 * N_pulse)
ru_vals = np.array(
    [autocorr_lag(best_pulse, int(round(tau / dt))) for tau in tau_range]
)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(t * 1e6, best_pulse / np.max(np.abs(best_pulse)), "-b", label="u(t) (norm)")
ax.plot(
    tau_range * 1e6,
    np.real(ru_vals) / np.max(np.abs(np.real(ru_vals))),
    "-r",
    label=r"$r_u(\tau)$ (norm)",
)
ax.set_xlabel("Time (us)")
ax.legend()
ax.grid(True)
ax.set_title("Pulse and autocorrelation")
fig.savefig(os.path.join(output_dir, "pulse_autocorrelation.png"), dpi=150)
plt.show()

# %%
# Plot: spectrum

U_f = np.fft.fftshift(np.fft.fft(best_pulse, n=N_fft)) * dt

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(f / 1e3, np.abs(U_f), "-b", label=r"$|\tilde{u}(f)|$")
ax.plot(f / 1e3, np.abs(U_f) ** 2, "-r", label=r"$|\tilde{u}(f)|^2$")
ax.set_xlabel("Frequency (kHz)")
ax.legend()
ax.grid(True)
ax.set_xlim([-1000, 1000])
ax.set_title("Pulse spectrum")
fig.savefig(os.path.join(output_dir, "pulse_spectrum.png"), dpi=150)
plt.show()

# %%
# Plot: spectrum vs mask

Upow_norm = np.abs(U_f) ** 2 / np.max(np.abs(U_f) ** 2)
sdB = 10 * np.log10(Upow_norm + 1e-30)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(f / 1e3, sdB, "-b", label=r"$10\log_{10}|\tilde{u}(f)|^2$")
ax.plot(freqm_kHz, mask_dB, "r-o", label="GSM Mask")
ax.set_xlabel("Frequency (kHz)")
ax.set_ylabel("dB")
ax.legend()
ax.grid(True)
ax.set_xlim([-1000, 1000])
ax.set_ylim([-80, 5])
ax.set_title("Spectrum vs Mask")
fig.savefig(os.path.join(output_dir, "spectrum_vs_mask.png"), dpi=150)
plt.show()

# %%
# Transmitted signal (Part 2.3)


# from digicomm_book.ipynb
def linearStreamModulator(ut, dk, numSamplePerSymbol):
    idx = np.arange(len(ut))
    xt = np.zeros(len(ut) + (len(dk) - 1) * numSamplePerSymbol, dtype=complex)
    for i_symbol in range(len(dk)):
        xt[idx + i_symbol * numSamplePerSymbol] += ut * dk[i_symbol]
    return xt


x_vec = np.array([1, -1, 1, 1, -1, 1, -1, -1, -1, 1])
nsps = int(round(best_tb / dt))

xt = linearStreamModulator(best_pulse, x_vec, nsps)
t_xt = np.arange(len(xt)) * dt * 1e6

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(t_xt, np.real(xt), "-b", label=r"$\Re[x(t)]$")
ax.set_xlabel("Time (us)")
ax.set_ylabel("Amplitude")
ax.legend()
ax.grid(True)
ax.set_title("Transmitted signal x(t)")
fig.savefig(os.path.join(output_dir, "transmitted_signal.png"), dpi=150)
plt.show()

# %%
# Matched filter output (Part 2.3)

mf = np.conj(best_pulse[::-1])
mf_out = np.convolve(xt, mf, mode="full") * dt
t_mf = np.arange(len(mf_out)) * dt * 1e6

# sample at symbol times
center = len(best_pulse) - 1
samp_idx = center + np.arange(len(x_vec)) * nsps
samp_idx = samp_idx[samp_idx < len(mf_out)]

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(t_mf, np.real(mf_out), "-b", label="MF output")
ax.plot(
    samp_idx * dt * 1e6,
    np.real(mf_out[samp_idx]),
    "ro",
    markersize=8,
    label="Sample points",
)
ax.set_xlabel("Time (us)")
ax.set_ylabel("Amplitude")
ax.legend()
ax.grid(True)
ax.set_title("Matched filter output")
fig.savefig(os.path.join(output_dir, "matched_filter_output.png"), dpi=150)
plt.show()

# %%
