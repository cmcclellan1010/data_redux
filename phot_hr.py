import numpy as np
from copy import deepcopy
from astropy.stats import mad_std
from astropy.io import fits
from photutils import DAOStarFinder
from photutils import aperture_photometry, CircularAperture
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

cwd = "/home/connor/programming/PycharmProjects/AST4723C/observing_project/hr/"
filter1 = 'v'
filter2 = 'i'


def load_fits(filename):
    hdu_list = fits.open(filename)
    data = hdu_list[0].data
    hdu_list.close()
    return data


def intersect_stars(star_table1, star_table2, da2):
    """Returns information about the stars in the star table closest to the specified (x, y) positions."""
    id_list = []
    coordlist = []
    for i in range(len(star_table1)):
        coordlist.append((np.float64(star_table1[i][1]), np.float64(star_table1[i][2])))

    for j in range(len(coordlist)):
        x = coordlist[j][0]
        y = coordlist[j][1]
        d_list = []
        for i in range(np.shape(star_table2)[0]):
            x2 = np.float64(star_table2[i][1])
            y2 = np.float64(star_table2[i][2])
            d = np.sqrt((x-x2)**2 + (y-y2)**2)
            d_list.append(d)
        id = d_list.index(np.min(d_list))
        dist = np.min(d_list)
        if dist < 0.5:
            id_list.append(id)
    trimmed = star_table2[id_list]
    centroids = []
    for i in range(len(trimmed)):
        centroids.append(((np.float64(trimmed[i][1]), np.float64(trimmed[i][2]))))
    source_aperture = CircularAperture(centroids, r=3.)
    source_table = aperture_photometry(da2, source_aperture)
    return source_table


def star_find(data, fwhm=12., threshold=50.):
    bkg_sigma = mad_std(data)
    if bkg_sigma == 0.0:
        thresh = 1000
    else:
        thresh = threshold*bkg_sigma
    daofind = DAOStarFinder(fwhm=fwhm, threshold=thresh)
    sources = daofind(data)
    return sources


def ap_phot(sources, data, source_r=3.):
    global fig
    centroids = (sources['xcentroid'], sources['ycentroid'])
    source_aperture = CircularAperture(centroids, r=source_r)
    source_area = source_aperture.area()
    source_table = aperture_photometry(data, source_aperture)

    sky_aperture = CircularAperture((290, 130), r=3.)
    sky_area = sky_aperture.area()
    sky_table = aperture_photometry(data, sky_aperture)

    sky_subtracted_source = deepcopy(source_table)

    for i in range(np.shape(centroids)[1]):
        sky_value = sky_table[0][3]
        sky_per_pix = sky_value / sky_area
        sky_behind_source = sky_per_pix * source_area
        sky_subtracted_source[i][3] -= sky_behind_source

    return sky_subtracted_source


def plot_phot(star_table, data, color, alpha=1):
    global fig, red_giant_starid, ref_starid, red_giant_iter, ref_star_iter, rg2_starid, rg2_star_iter
    centroids = (star_table['xcenter'], star_table['ycenter'])
    source_aperture = CircularAperture(centroids, r=3)

    sky_aperture = CircularAperture((290, 130), r=3.)
    sky_table = aperture_photometry(data, sky_aperture)

    plt.imshow(data, cmap='gray', origin='lower', vmin=0, vmax=1500, alpha=alpha)
    for i in range(np.shape(centroids)[1]):
        if (np.float64(centroids[0][i]) - 739.5)**2 + (np.float64(centroids[1][i]) - 699.)**2 < 0.5:
            plt.annotate("NGC 7078 682",
                         xy=(np.float64(star_table[i][1]) + 5.,
                             np.float64(star_table[i][2]) + 5.),
                         color=color, alpha=alpha)
            red_giant_starid = int(star_table[i][0])
            red_giant_iter = i
        elif (np.float64(centroids[0][i]) - 859.)**2 + (np.float64(centroids[1][i]) - 788.)**2 < 0.5:
            plt.annotate("NGC 7978 341",
                         xy=(np.float64(star_table[i][1]) + 5.,
                             np.float64(star_table[i][2]) + 5.),
                         color=color, alpha=alpha)
            ref_starid = int(star_table[i][0])
            ref_star_iter = i
        elif (np.float64(centroids[0][i]) - 500.) ** 2 + (np.float64(centroids[1][i]) - 740.3) ** 2 < 0.5:
            plt.annotate("NGC 7978 1030",
                         xy=(np.float64(star_table[i][1]) + 5.,
                             np.float64(star_table[i][2]) + 5.),
                         color=color, alpha=alpha)
            rg2_starid = int(star_table[i][0])
            rg2_star_iter = i
        else:
            plt.annotate(str(star_table[i][0]),
                         xy=(np.float64(star_table[i][1]) + 5.,
                             np.float64(star_table[i][2]) + 5.),
                         color=color, alpha=alpha)
    source_aperture.plot(color=color, lw=1.5, alpha=0.5)
    plt.annotate("SKY",
                 xy=(np.float64(sky_table[0][1]) + 15.,
                     np.float64(sky_table[0][2]) + 15.),
                 color="white")
    sky_aperture.plot(color="orange", lw=0.5, alpha=0.5)
    plt.tight_layout()


def column(matrix, i, dtype=np.float64):
    return [dtype(row[i]) for row in matrix]


def find_ref(data):
    star_ids = column(data, 0, dtype=np.int)
    for i in range(len(data)):
        if (column(data, 1)[i] - 859) ** 2 + (column(data, 2)[i] - 788) ** 2 < 0.5:
            id = star_ids.index(data[i][0])
    return data[id]


# Make the center mask
h = 454.
k = 778.
r = 73.


def maskf(x, y, (y1, y2), (x1, x2)):
    return x1 < x < x2 and y1 < y < y2


def check_if_mask(x, y):
    if (x - h)**2 + (y - k)**2 < r**2:
        mask = True
    elif maskf(x, y, (853, 903), (403, 415)):
        mask = True
    elif maskf(x, y, (725, 736), (366, 376)):
        mask = True
    elif maskf(x, y, (688, 885), (464, 472)):
        mask = True
    elif maskf(x, y, (785, 798), (523, 530)):
        mask = True
    elif maskf(x, y, (701, 709), (433, 439)):
        mask = True
    elif maskf(x, y, (614, 628), (561, 574)):
        mask = True
    elif maskf(x, y, (600, 632), (340, 351)):
        mask = True
    elif maskf(x, y, (804, 814), (340, 348)):
        mask = True
    elif maskf(x, y, (714, 724), (375, 383)):
        mask = True
    else:
        mask = False
    return mask


da1 = fits.open('/home/connor/programming/PycharmProjects/AST4723C/observing_project/hr/m15_v_30s.fit')[0].data
da2 = fits.open('/home/connor/programming/PycharmProjects/AST4723C/observing_project/hr/m15_i_30s.fit')[0].data

# mask the center region
dimensions = np.shape(da1)
mask = np.ones(dimensions)

for i in range(dimensions[0]):
    for j in range(dimensions[1]):
        if check_if_mask(float(i), float(j)):
            mask[i][j] = 0
            da1[i][j] = 0.0
            da2[i][j] = 0.0

sources1 = star_find(da1)
star_table1 = ap_phot(sources1, da1)
fig = plt.figure(figsize=(17, 17))
plot_phot(star_table1, da1, 'blue')
plt.show()
print "Stars found in image 1: ", len(star_table1)

sources2 = star_find(da2)
star_table2 = ap_phot(sources2, da2)
fig = plt.figure(figsize=(17, 17))
plot_phot(star_table2, da2, 'red')
plt.show()
print "Stars found in image 2: ", len(star_table2)

phot_2 = intersect_stars(star_table1, star_table2, da2)
print "Stars in both images: ", len(phot_2)

phot_1 = intersect_stars(star_table2, star_table1, da1)
print "Stars in both images: ", len(phot_1)


fig = plt.figure(figsize=(17, 17))
plot_phot(phot_1, da1, color="blue")
plot_phot(phot_2, da2, color="red", alpha=0.5)
plt.show()


def flux(counts):
    return np.float64(counts)*2.54/30.      # in electrons per second


def flux_err(counts):
    return flux(np.sqrt(counts))


Vref_flux = flux(find_ref(phot_1)[3])
Vref_flux_err = flux_err(find_ref(phot_1)[3])
print Vref_flux, Vref_flux_err

Iref_flux = flux(find_ref(phot_2)[3])
Iref_flux_err = flux_err(find_ref(phot_2)[3])
print Iref_flux, Iref_flux_err

VZ = 12.54 + 2.5*np.log10(Vref_flux)
VZ_err = np.sqrt(0.01**2 + ((2.5/np.log(10))*Vref_flux_err/Vref_flux)**2)
print VZ, VZ_err

VI = 11.22 + 2.5*np.log10(Iref_flux)
VI_err = np.sqrt(0.01**2 + ((2.5/np.log(10))*Iref_flux_err/Iref_flux)**2)
print VI, VI_err


def flux_to_mag(aflux, zpt):
    return -2.5*np.log10(aflux) + zpt


def flux_to_mag_err(aflux, aflux_err, zpt_err):
    return np.sqrt(((2.5/np.log(10))*aflux_err/aflux)**2 + zpt_err**2)


mags = []
mags_err = []
v_minus_i = []
v_minus_i_err = []
for i in range(len(phot_1)):
    v_flux = flux(phot_1[i][3])
    v_flux_err = flux_err(phot_1[i][3])
    i_flux = flux(phot_2[i][3])
    i_flux_err = flux_err(phot_2[i][3])

    mags.append(flux_to_mag(v_flux, VZ))
    mags_err.append(flux_to_mag_err(v_flux, v_flux_err, VZ_err))

    v_minus_i.append(flux_to_mag(v_flux, VZ) - flux_to_mag(i_flux, VI))
    v_minus_i_err.append(np.sqrt(flux_to_mag_err(v_flux, v_flux_err, VZ_err)**2 + flux_to_mag_err(i_flux, i_flux_err, VI_err)**2))


def convert_to_absolute(relative):
    d = 10000.
    mu = 5.*np.log10(d) - 5.
    absolute = relative - mu
    return absolute


v_mag, v_mag_err = np.zeros(len(phot_1)), np.zeros(len(phot_1))
vi_mag, vi_mag_err = np.zeros(len(phot_1)), np.zeros(len(phot_1))


for i in range(len(phot_1)):
    v_mag[i] = convert_to_absolute(mags[i])
    v_mag_err[i] = mags_err[i]
    vi_mag[i] = v_minus_i[i]
    vi_mag_err[i] = v_minus_i_err[i]
    print phot_1[i][0], '{0:.4f}'.format(v_mag[i]), '{0:.4f}'.format(v_mag_err[i]), '{0:.4f}'.format(vi_mag[i]), '{0:.4f}'.format(vi_mag_err[i])

plt.figure(figsize=(17, 17))
plt.errorbar(vi_mag, v_mag, xerr=vi_mag_err, yerr=v_mag_err, fmt='ko', markersize=3, ecolor='k', elinewidth=1, capsize=2, alpha=0.5)


for i in range(len(phot_1)):
    if phot_1[i][0] == red_giant_iter:
        plt.annotate("NGC 7078 682",
                     xy=(vi_mag[i]+.0001, v_mag[i]+.0001),
                     color='k',
                     fontsize=8,
                     fontweight='bold')
    elif phot_1[i][0] == ref_star_iter:
        plt.annotate("NGC 7978 341",
                     xy=(vi_mag[i]+.0001, v_mag[i]+.0001),
                     color='k',
                     fontsize=8,
                     fontweight='bold')
    elif phot_1[i][0] == rg2_star_iter:
        plt.annotate("NGC 7978 1030",
                     xy=(vi_mag[i]+.0001, v_mag[i]+.0001),
                     color='k',
                     fontsize=8,
                     fontweight='bold')
    else:
        plt.annotate(phot_1[i][0],
                     xy=(vi_mag[i]+.0001, v_mag[i]+.0001),
                     color='k',
                     fontsize=8,
                     fontweight='bold')


plt.xlabel("V - I (mag)")
plt.ylabel("V (abs. mag)")
plt.gca().invert_yaxis()
plt.show()

print red_giant_starid
print red_giant_iter


