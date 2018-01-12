import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import itertools
from copy import deepcopy


def rms(x):
    return (np.absolute(np.mean(x**2) - (np.mean(x))**2))**0.5


def shuffle_in_unison(a, b):
    rng_state = np.random.get_state()
    np.random.shuffle(a)
    np.random.set_state(rng_state)
    np.random.shuffle(b)


filter = 'r'
filename = "/home/connor/programming/PycharmProjects/AST4723C/observing_project/variability/photometry_"+filter+".fits"
filename2 = "/home/connor/programming/PycharmProjects/AST4723C/observing_project/variability/photometry_"+filter+"_CD.fits"
da = fits.open(filename)[1].data
da2 = fits.open(filename2)[1].data

# convert counts to e-/s
gain = 1.49         # electrons per ADU
exposure_time = 30  # seconds

flux1_list = []
fluxA_list = []
fluxB_list = []
fluxC_list = []
fluxD_list = []
mjd_list = []

for i in range(len(da)):
    flux1_list.append(da[i][1]*gain/exposure_time)
    fluxA_list.append(da[i][2]*gain/exposure_time)
    fluxB_list.append(da[i][3]*gain/exposure_time)
    fluxC_list.append(da2[i][2]*gain/exposure_time)
    fluxD_list.append(da2[i][3]*gain/exposure_time)
    mjd_list.append(da[i][4])

target = np.asarray(flux1_list)
refA = np.asarray(fluxA_list)
refB = np.asarray(fluxB_list)
refC = np.asarray(fluxC_list)
refD = np.asarray(fluxD_list)
mjd = np.asarray(mjd_list)

# sort the data by MJD
lists = sorted(itertools.izip(*[mjd, target, refA, refB, refC, refD]))
mjd, target, refA, refB, refC, refD = np.asarray(list(itertools.izip(*lists)))
hr = (mjd - mjd[0])*24.
hra = range(len(hr))

r_avg = []
r_avg_err = []
for i in range(len(hr)):
    ref_fluxes = [refA[i], refC[i], refD[i]]
    r_avg.append(np.mean(ref_fluxes))
r_avg = np.array(r_avg)


def get_residuals(refstar, avg):
    reslist = []
    for i in range(len(refstar)):
        div = refstar[i]/avg[i]
        mean = np.mean(refstar/avg)
        res = div - mean
        reslist.append(res)
    return np.array(reslist)


A_res = get_residuals(refA, r_avg)
C_res = get_residuals(refC, r_avg)
D_res = get_residuals(refD, r_avg)

for i in range(len(A_res)):
    r_avg_err.append(rms(np.array([A_res[i], C_res[i], D_res[i]]))/np.sqrt(3))
r_avg_err = np.array(r_avg_err)

t_div = target/r_avg
t_err = np.sqrt(target)

# find the target uncertainty
t_div_err = t_div*np.sqrt((t_err/target)**2 + (r_avg_err/r_avg)**2)

from astropy.stats import LombScargle
radius = 0.01
upper_range = 0.054105527 + radius
lower_range = 0.054105527 - radius

frequency = np.linspace(lower_range, upper_range, 1000)
power = LombScargle(hr, t_div, t_div_err).power(frequency)

print "Power at peak frequency: ", np.max(power)
print "Peak frequency: ", frequency[power.argmax()]

# PLOT LOMB-SCARGLE
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12.8, 7.6))
ax1.plot(frequency, power, 'k')
ax1.set_title("Lomb-Scargle Periodogram of HD 345439")
ax1.set_xlabel("Frequency (Cycles per hour)")
ax1.set_ylabel("Power")
ax1.set_xlim([np.min(frequency), np.max(frequency)])
ax1.annotate("Power at peak frequency: %.4f" % np.max(power), xy=(0.01, 0.85), xycoords='axes fraction')
ax1.annotate("Peak Frequency: %.11f" % frequency[power.argmax()], xy=(0.01, 0.90), xycoords='axes fraction')

period = 1./frequency[power.argmax()]

period_in_s = 60.*60.*period
m, s = divmod(period_in_s, 60.)
h, m = divmod(m, 60.)
ax1.annotate("Period: %.02d h %.02d m %.5f s" % (h, m, s), xy=(0.01, 0.75), xycoords='axes fraction')

best_frequency = frequency[np.argmax(power)]
hr_fit = np.linspace(0, period)
t_div_fit = LombScargle(hr, t_div, t_div_err).model(hr_fit, best_frequency)
ax2.plot(hr_fit/period, t_div_fit, 'limegreen')

phased_hr = []
for i in range(len(hr)):
    phased_hr.append(divmod(hr[i], period)[1]/period)
phased_hr = np.array(phased_hr)
ax2.errorbar(phased_hr, t_div, yerr=t_div_err, fmt='ko', markersize=3, ecolor='k', elinewidth=1, capsize=2, alpha=0.5)
ax2.set_title("Phased data at frequency = %.6f" % best_frequency)
ax2.set_ylabel("Target / R average Flux")
ax2.set_xlabel("Phase")
ax2.set_xlim([np.min(phased_hr)-0.1, np.max(phased_hr)+0.1])
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.97, wspace=.2, hspace=.25)
plt.show()

# MONTE-CARLO UNCERTAINTY
n_tries = 10000
powers = []
similars = []

hits = 0
for i in range(n_tries):
    # scramble the data
    t_scramble = deepcopy(t_div)
    err_scramble = deepcopy(t_div_err)

    shuffle_in_unison(t_scramble, err_scramble)

    mfrequency = np.linspace(lower_range, upper_range, 1000)
    mpower = LombScargle(hr, t_scramble, err_scramble).power(frequency)
    powers.append(np.max(mpower))

    if np.max(mpower) > np.max(power):
        hits += 1

print "False detections: ", hits
print "Percent false: ", 100*float(hits)/float(len(powers))
print "Percent confidence: ", 100. - 100*float(hits)/float(len(powers))