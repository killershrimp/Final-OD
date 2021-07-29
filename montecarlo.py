from libs.odlib import *
import matplotlib.pyplot as plt
import numpy as np

data = parseMonteCarloData("obs.txt")


k_iter = 100_000
orb_els = []
k_des_JD = timeToJD("03:18:49") + getJD(2021, 7, 11)    # middle observation

samples = 0
for iteration in range(k_iter):
    # print(iteration, "before")
    obs = []        # not including MA
    mas = []
    for i in range(3):
        date = data[i][0]
        ra = data[i][1]
        dec = data[i][2]
        ra_unc = data[i][3]
        dec_unc = data[i][4]
        samp_ra = np.random.normal(ra, ra_unc)
        samp_dec = np.random.normal(dec, dec_unc)
        # samp_ra = ra
        # samp_dec = dec
        e_s_vec = data[i][5]
        obs.append((date, samp_ra, samp_dec, e_s_vec))

    # o_obs = parseODData("YunMoGInp.txt")
    # o_mog_ret = MOG(o_obs, log=False, want_ecl=True)
    # if o_mog_ret is None:
    #     print("Didn't converge!")
    #     exit(1)
    # o_pos_vec = o_mog_ret[0]
    # o_vel_vec = o_mog_ret[1]
    # o_orb_el_JD = o_mog_ret[2]
    #
    # o_des_JD = getJD(2021, 7, 17) + timeToJD("03:00:13")
    # o_orb_el_JD = incrementMA(o_orb_el_JD, o_des_JD)


    mog_ret = MOG(obs, want_ecl=True)
    if mog_ret is None:     # didn't converge
        # print("didn't converge")
        continue
    samples += 1
    r_vec = mog_ret[0]
    r2dot_vec = mog_ret[1]
    orb_el = mog_ret[2]

    orb_el = incrementMA(orb_el, k_des_JD)
    orb_els.append([orb_el[0].a, orb_el[0].e, math.degrees(orb_el[0].i), math.degrees(orb_el[0].APE), math.degrees(orb_el[0].LOAN), math.degrees(orb_el[0].MA)])

    # print(iteration, "after")

orb_els = np.array(orb_els)
k_bins = 100
rows = 3
cols = 2
# fig, axs = plt.subplots(rows, cols, figsize=(8, 6), tight_layout=True)

jpl_orb_el = [2.332866366174863E+00, 4.580825245891841E-01, 5.721647546127751E+00, 5.028093353213583E+01, 3.088649396519276E+02, 3.381874492423444E+02]


titles = ["Semi-major Axis", "Eccentricity", "Inclination", "Argument of Perihelion",
          "Longitude of Ascending Node", "Mean Anomaly at Middle Observation"]
x_axes = ["a (AU)", "e", "i (deg)", "ω (deg)", "Ω (deg)", "Mean Anomaly (deg)"]
# fig.canvas.draw()
for counter in range(6):
    plt.title(titles[counter] + " Distributions for " + str(samples) + " Samples:")
    plt.xlabel(x_axes[counter])
    plt.ylabel("Frequency")
    el_data = orb_els[:,counter]
    n, bins, patches = plt.hist(el_data, bins=k_bins, facecolor='#2ab0ff', edgecolor='#e0e0e0', linewidth=0.5, alpha=0.7)
    n = n.astype("int")
    for z in range(len(patches)):
        patches[z].set_facecolor(plt.cm.viridis(n[z] / max(n)))

    mean = np.mean(el_data)
    print(titles[counter])
    print("\t", mean)
    std = np.std(el_data)
    plt.axvline(mean, color="k", linestyle="dashed", linewidth=1, label="Mean: {:.5f}".format(mean))
    plt.axvline(mean-std, color="k", linewidth=1, label="Std. Dev.: {:.5f}".format(std))
    plt.axvline(mean+std, color="k", linewidth=1)
    plt.axvline(mean-2*std, color="k", linewidth=1)
    plt.axvline(mean+2*std, color="k", linewidth=1)
    plt.axvline(jpl_orb_el[counter], color="r", linewidth=1, label="JPL Estimate: {:.5f}".format(jpl_orb_el[counter]))

    plt.legend()
    plt.show()
plt.show()

