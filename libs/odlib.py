from libs.orbital_elements import OrbitalElements
import math
import numpy as np


# CONSTANTS
k_sol_days_per_yr = 365.2563835
k_gauss_days_per_yr = 2 * math.pi
k_gauss_days_per_sol = k_gauss_days_per_yr / k_sol_days_per_yr
k_km_per_AU = 1.496 * 1e8

k_grav = 6.67408 * 1e-11        # m^3 / kg / s^2
k_c = 2.99792 * 1e5 * 1e3       # km / s
k_c_AU = 173.144643267          # AU / mean solar day
k_m_sun = 1.9891 * 1e30         # kg
k_m_earth = 5.972 * 1e24        # kg

k_obliquity = math.radians(23.4374) # Earth obliquity

k_k = 2 * math.pi / k_sol_days_per_yr   # Gauss's constant
# k_k = 0.0172020989484       # Gaussian constant
k_mu = k_k ** 2


# METHODS
def inRange(a, b, range):
    return abs(a - b) <= range


def readVectorsData(filename):
    """
    :param filename: file name
    :return: big tuple of tuple (JD, np.array pos, np.array vel)
    """
    ret = []
    with open(filename, "r") as file:
        line = file.readline().strip()
        while line != "":
            JD_cutoff = line.find("=")
            JD = float(line[:JD_cutoff].strip())
            raw_pos = file.readline().strip()
            raw_vel = file.readline().strip()
            raw_lt_rg_rr = file.readline().strip()  # not used to test

            pos = []
            x_str_pos = raw_pos.find("X")
            y_str_pos = raw_pos.find("Y")
            z_str_pos = raw_pos.find("Z")

            pos.append(float(raw_pos[x_str_pos + 3:y_str_pos].strip()))  # x comp
            pos.append(float(raw_pos[y_str_pos + 3:z_str_pos].strip()))  # y comp
            pos.append(float(raw_pos[z_str_pos + 3:].strip()))  # z comp

            vel = []
            dx_str_pos = raw_vel.find("VX")
            dy_str_pos = raw_vel.find("VY")
            dz_str_pos = raw_vel.find("VZ")

            vel.append(float(raw_vel[dx_str_pos + 3:dy_str_pos].strip()) / k_gauss_days_per_sol)
            vel.append(float(raw_vel[dy_str_pos + 3:dz_str_pos].strip()) / k_gauss_days_per_sol)
            vel.append(float(raw_vel[dz_str_pos + 3:].strip()) / k_gauss_days_per_sol)

            ret.append((JD, np.array(pos), np.array(vel)))
            line = file.readline().strip()
    return ret


def readOrbElData(filename):
    """
    :param filename: name of file containing orbital elements
    :return: big tuple of tuple (a, e, i, LOAN, TA), with angles in rad
    """
    ret = []
    with open(filename, "r") as file:
        line = file.readline().strip()
        while line != "":
            JD_cutoff = line.find("=")
            JD = float(line[:JD_cutoff].strip())
            raw_EC_QR_IN = file.readline().strip()
            raw_OM_W_TP = file.readline().strip()
            raw_N_MA_TA = file.readline().strip()
            raw_A_AD_PR = file.readline().strip()

            e_str_pos = raw_EC_QR_IN.find("EC")
            qr_str_pos = raw_EC_QR_IN.find("QR")
            i_str_pos = raw_EC_QR_IN.find("IN")

            e = float(raw_EC_QR_IN[e_str_pos + 3:qr_str_pos].strip())
            i = float(raw_EC_QR_IN[i_str_pos + 3:].strip())

            LOAN_str_pos = raw_OM_W_TP.find("OM")
            w_str_pos = raw_OM_W_TP.find("W")
            TP_str_pos = raw_OM_W_TP.find("Tp")

            LOAN = float(raw_OM_W_TP[LOAN_str_pos + 3:w_str_pos].strip())
            APE = float(raw_OM_W_TP[w_str_pos + 3:TP_str_pos].strip())
            PET = float(raw_OM_W_TP[TP_str_pos + 3:].strip())

            TA_str_pos = raw_N_MA_TA.find("TA")
            # TA = float(raw_N_MA_TA[TA_str_pos + 3:].strip())
            MA_str_pos = raw_N_MA_TA.find("MA")
            MA = float(raw_N_MA_TA[MA_str_pos + 3:TA_str_pos])

            a_str_pos = raw_A_AD_PR.find("A")
            AD_str_pos = raw_A_AD_PR.find("AD")
            a = float(raw_A_AD_PR[a_str_pos + 3:AD_str_pos].strip())

            ret.append((JD, OrbitalElements(a, e, math.radians(i), math.radians(LOAN), math.radians(APE), math.radians(MA), PET)))
            line = file.readline().strip()
    return ret


def getCorrectQuadrant(sin_t, cos_t):
    """
    Given sine and cosine of angle t, return t bounded between 0 and 2pi
    """
    if sin_t >= 0 and cos_t >= 0:
        # 1st quadrant
        return math.asin(sin_t)
    if sin_t >= 0 and cos_t < 0:
        # 2nd quadrant
        return math.acos(cos_t)
    if cos_t >= 0:
        # 4th quadrant
        return 2*math.pi + math.asin(sin_t)
    # 3rd quadrant
    return math.pi - math.asin(sin_t)


def bound0To2Pi(angle):
    """
    :param angle: angle in radians
    :return: angle bound between 0 and 2pi rads
    """
    while angle < 0:
        angle += 2 * math.pi
    return angle % (2 * math.pi)


def parseODData(filename):
    """
    File format: for each night
        UT date/time at middle of observation (JD)
        RA (sex. hr.), DEC (sex. deg)
        Earth->Sun vector (AU, eq cartesian, J2000, apparent)
    :param filename: file name
    :return: (obs date/time, (RA, DEC), e_s_vec)
    """
    obs = []
    with open(filename, "r") as file:
        line = file.readline()
        while line != "":
            line = line.split(" ")
            date = getJD(int(line[0]), int(line[1]), int(line[2])) + timeToJD(line[3])
            ra_sex = hrsToDecDeg(line[4])
            dec_sex = sexDegToDecDeg(line[5])
            earth_sun_vec = np.array([
                float(line[6]),
                float(line[7]),
                float(line[8])
            ])
            obs.append((date, ra_sex, dec_sex, earth_sun_vec))
            line = file.readline()
    return obs


def timeToJD(time):
    """
    Converts time string to JD
    :param time: "HH:MM:SS.S"
    :return: decimal day
    """
    hrs, min, sec = time.split(":")
    tot_sec = int(hrs) * 3600 + int(min) * 60 + float(sec)
    return tot_sec / (24 * 3600)


def parseMonteCarloData(filename):
    """
    File format: for each night
        UT date/time at middle of observation (JD)
        RA (sex. hr.), DEC (sex. deg)
        RA uncert., Dec. uncert.
        Earth->Sun vector (AU, eq cartesian, J2000, apparent)
    :param filename: file name
    :return: (obs date/time, (RA, DEC), e_s_vec)
    """
    obs = []
    with open(filename, "r") as file:
        line = file.readline()
        while line != "":
            line = line.split(" ")
            date = getJD(int(line[0]), int(line[1]), float(line[2])) + timeToJD(line[3])
            ra_sex = hrsToDecDeg(line[4])
            dec_sex = sexDegToDecDeg(line[5])
            ra_uncert = float(line[6])
            dec_uncert = float(line[7])
            earth_sun_vec = np.array([
                float(line[8]),
                float(line[9]),
                float(line[10])
            ])
            obs.append((date, ra_sex, dec_sex, ra_uncert, dec_uncert, earth_sun_vec))
            line = file.readline()
    return obs

def getJD(Y, M, D):
    j0 = 367*Y - int(7 * (Y + int((M+9)/12)) / 4) + int(275*M/9) + D + 1721013.5
    return j0


def decDegToSexDeg(deg, round_places=None):
    """
    Converts angles from decimal degrees to sexigesimal degrees
    :param deg: decimal degrees
    :param round_places: Either
        - None: no rounding will be done
        - [a number]: answer will be rounded to that decimal place
    :return: angle in dec deg as a string (digits separated by ':')
    """
    remaining = abs(deg)
    hrs = int(remaining)
    remaining -= hrs
    min = int(remaining * 60)
    remaining -= min / 60
    sec = remaining * 3600
    if round_places is not None:
        sec = round(sec, round_places)
    if hrs < 10:
        hrs = "0" + str(hrs)
    if min < 10:
        min = "0" + str(min)
    if sec < 10:
        sec = "0" + str(sec)
    return ("-" if deg < 0 else "+") + str(hrs) + ":" + str(min) + ":" + str(sec)


def decDegToHr(deg, round_places=None):
    """
    Converts angles from decimal degrees to hrs
    :param deg: decimal degrees
    :param round_places: Either
        - None: no rounding will be done
        - [a number]: answer will be rounded to that decimal place
    :return: angle in hrs as a string (digits separated by ':')
    """
    deg *= 24/360
    hrs = int(deg)
    deg -= hrs
    min = int(deg * 60)
    deg -= min / 60
    sec = deg * 3600
    if round_places is not None:
        sec = round(sec, round_places)
    if hrs < 10:
        hrs = "0" + str(hrs)
    if min < 10:
        min = "0" + str(min)
    if sec < 10:
        sec = "0" + str(sec)
    return str(hrs) + ":" + str(min) + ":" + str(sec)


def hrsToDecDeg(coord):
    """
    :param coord: string of sexagesimal angle measure in hrs, separated by ":"
    :return: angle in decimal degrees
    """
    chars = coord.split(":")
    sign = np.sign(float(chars[0]))
    hrs = abs(float(chars[0])) * 3600
    min = abs(float(chars[1])) * 60
    sec = abs(float(chars[2]))
    arcsecs = 360/24 * (hrs + min + sec)
    return sign * arcsecToDecDeg(arcsecs)


def sexDegToDecDeg(coord):
    """
    :param coord: string of sexagesimal angle measure in deg, separated by ":"
    :return: angle in decimal deg
    """
    chars = coord.split(":")
    arcsecs = np.sign(float(chars[0])) * ((abs(float(chars[0])) * 60 + float(chars[1])) * 60 + float(chars[2]))
    return arcsecToDecDeg(arcsecs)


def arcsecToDecDeg(num):
    return num / 3600


def incrementMA(orb_el_JD, des_JD):
    """
    Get Mean Anomaly for certain day based on current MA
    :param orb_el_JD: (OrbitalElements object, decimal JD)
    :param des_JD: desired JD (as decimal)
    :returns:
        - OrbitalElements obj with adjusted Mean Anomaly
        - Desired JD
    """
    old_JD = orb_el_JD[1]
    orb_el = orb_el_JD[0]

    n = orb_el.a ** (-3 / 2) * k_k
    d_JD = des_JD - old_JD
    orb_el.MA += n * d_JD
    return orb_el, des_JD


def getHVec(pos, vel):
    """
    :param pos: pos vec (AU)
    :param vel: vel vec (AU/day)
    :return: angular momentum vector rounded to 6 decimal places
    """
    prod = np.cross(pos, vel)
    ret = (round(prod[0], 6), round(prod[1], 6), round(prod[2], 6))
    return np.array(ret)


def getA(r_vec, v_vec):
    """
    :param r_vec: pos vec (AU)
    :param v_vec: vel vec (AU/day)
    :return: semi-major axis in AU
    """
    r_mag = np.linalg.norm(r_vec)
    v_mag = np.linalg.norm(v_vec)
    return 1 / (2 / r_mag - v_mag ** 2)


def getE(h_vec, a):
    """
    :param h_vec: ang momentum vector
    :param a: semi-major axis (AU)
    :return: eccentricity
    """
    h_mag = np.linalg.norm(h_vec)
    return (1 - h_mag**2 / a) ** (1/2)


def getI(h_vec):
    """
    :param h_vec: ang momentum vector
    :return: inclination (rad)
    """
    hx = h_vec [0]
    hy = h_vec[1]
    hz = h_vec[2]
    return math.atan2((hx ** 2 + hy ** 2) ** (1/2), hz)


def getU(r_vec, i, LOAN):
    """
    :param r_vec: pos vec (AU)
    :param i: inclination (rad)
    :param LOAN: long. of asc. node (rad)
    :return: true longitude (rad)
    """
    x = r_vec[0]
    y = r_vec[1]
    z = r_vec[2]
    r_mag = np.linalg.norm(r_vec)
    sin_u = z / (r_mag * math.sin(i))
    cos_u = (x * math.cos(LOAN) + y * math.sin(LOAN)) / r_mag
    return getCorrectQuadrant(sin_u, cos_u)


def getAPE(r_vec, i, LOAN, TA):
    """
    :param r_vec: pos vec (AU)
    :param i: inclination (rad)
    :param LOAN: long. of asc. node (rad)
    :param TA: true anomaly (rad)
    :return: argument of perihelion (rad)
    """
    U = getU(r_vec, i, LOAN)
    return bound0To2Pi(U - TA)


# Long. of Asc. Node
def getLOAN(h_vec, i):  # sorry for reminding you about your student loans, TA :(
    """
    :param h_vec: ang momentum vec
    :param i: inclination (rad)
    :return: Longitude of Ascending Node (rad)
    """
    cos_Omega = - h_vec[1] / (math.sin(i) * np.linalg.norm(h_vec))
    sin_Omega = h_vec[0] / (math.sin(i) * np.linalg.norm(h_vec))
    return getCorrectQuadrant(sin_Omega, cos_Omega)


def getTA(a, e, h_vec, r_vec, v_vec):
    """
    :param a: semi-major axis (AU)
    :param e: eccentricity
    :param h_vec: ang momentum vec
    :param r_vec: pos vec (AU)
    :param v_vec: vel vec (AU/day)
    :return: true anomaly (rad)
    """
    h_mag = np.linalg.norm(h_vec)
    r_mag = np.linalg.norm(r_vec)
    sin_ta = a * (1 - e**2) / (e * h_mag) * (r_vec @ v_vec) / r_mag
    cos_ta = (a * (1 - e**2) / r_mag - 1) / e
    return getCorrectQuadrant(sin_ta, cos_ta)


def getEccAnom(r_vec, a, e, TA):
    """
    :param r_vec: position vector (AU)
    :param a: semi-major axis (AU)
    :param e: eccentricity
    :param TA: true anomaly (rad)
    :return: Eccentric Anomaly (rad)
    """
    r_mag = np.linalg.norm(r_vec)
    cos_E = (1 - r_mag / a) / e
    acos_E = math.acos(cos_E)
    if TA >= math.pi:
        return 2 * math.pi - acos_E
    return acos_E


def getMA(ecc_anom, e):
    """
    :param ecc_anom: eccentric(rad)
    :param e: eccentricity
    :return: mean anomaly (rad)
    """
    return ecc_anom - e * math.sin(ecc_anom)


def getPET(a, MA, JD):
    """
    :param a: semi-major axis (AU)
    :param MA: mean anomaly
    :param JD: current julian date
    :return: Time of Perihelion Passage (decimal JD)
    """
    denom = (k_k ** 2 / a ** 3) ** (1/2)
    return JD - MA / denom


def getInfinityStones(JD, r_vec, v_vec, want_tup=False):
    """
    Get 6 orbital elements (actually 7; gives both mean anomaly and last time of peri passage)
    :param r_vec: position vector
    :param v_vec: velocity vector
    :return: a, e, i, LOAN, APE, MA, PET, JD
    """
    # get orbital elements
    h_vec = getHVec(r_vec, v_vec)
    a = getA(r_vec, v_vec)
    e = getE(h_vec, a)
    i = getI(h_vec)
    LOAN = getLOAN(h_vec, i)
    TA = getTA(a, e, h_vec, r_vec, v_vec)
    APE = getAPE(r_vec, i, LOAN, TA)
    E = getEccAnom(r_vec, a, e, TA)
    MA = getMA(E, e)
    PET = getPET(a, MA, JD)
    if want_tup:
        return (a, e, i, LOAN, APE, MA, PET), JD
    return OrbitalElements(a, e, i, LOAN, APE, MA, PET), JD


def getInitialFG(tau, r2_vec):
    """
    Implement eq 104 to get fi, gi using 2nd degree Taylor Series
    :param tau: t_i
    :param r2_vec: r2 vector
    :return: f_deg and g_deg
    """
    r2_mag = np.linalg.norm(r2_vec)
    f_deg = 1 - 1 / (2 * r2_mag ** 3) * tau ** 2
    g_deg = tau - 1 / (6 * r2_mag ** 3) * tau ** 3
    return f_deg, g_deg


def getPossibleRho2s(r2s, A, B):
    """
    Implement eq. 128
    """
    rho2s = []
    for r2 in r2s:
        rho2 = A + B / (r2 ** 3)
        rho2s.append(rho2)
    return rho2s


def SEL(A, B, E, F):
    """
    Use eq. 129 to get np array of possible r2s
    """
    # coefficients of SEL
    a = - (A ** 2 + A * E + F)
    b = - (2 * A * B + B * E)       # no mu bc Gunits
    c = - B ** 2                    # no mu bc Gunits

    coeff = [1, 0, a, 0, 0, b, 0, 0, c]
    coeff.reverse()
    return np.polynomial.polynomial.polyroots(coeff)


def getInitR2sRho2s(taus, sun_pos_vec, rhohat2, d):
    """
    Implement eqs 129 to get initial r2
    :param taus: Gtime intervals (tau, tau1, tau3)
    :param sun_pos_vec: pos vec Sun -> asteroid
    :param rhohat2: unit vec Earth -> asteroid
    :param d: Coeff from scalar equations of range (D0, D21, D22, D23)
    :return: possible r2s, corresponding rhos
    """
    sun_dist_mag = np.linalg.norm(sun_pos_vec)
    tau = taus[0]
    tau1 = taus[1]
    tau3 = taus[2]
    A1 = tau3 / tau
    B1 = A1 / 6 * (tau ** 2 - tau3 ** 2)
    A3 = -tau1 / tau
    B3 = A3 / 6 * (tau ** 2 - tau1 ** 2)

    A = - (A1 * d[1] - d[2] + A3 * d[3]) / d[0]
    B = - (B1 * d[1] + B3 * d[3]) / d[0]
    E = -2 * (rhohat2 @ sun_pos_vec)
    F = sun_dist_mag ** 2

    roots = SEL(A, B, E, F)
    isreal_root = np.isreal(roots)
    real_roots = []
    for i in range(len(isreal_root)):
        if isreal_root[i] and roots[i] > 0:
            real_roots.append(roots[i].real)

    real_roots.sort()
    rhos = getPossibleRho2s(real_roots, A, B)

    final_roots = []
    final_rhos = []
    for i in range(len(real_roots)):
        if rhos[i] >= 0 and real_roots[i] >= 0:
            final_rhos.append(rhos[i])
            final_roots.append(real_roots[i])

    return final_roots, final_rhos


def getFGTaylor(tau, r_vec, r_vec_dot, want_third):
    """
    Implements eqs. 109, 110 to get fi and gi using Taylor Series
    :param tau: gaussian dt from tau2
    :param r_vec: r vector
    :param r_vec_dot: r dot vector
    :param want_third: (T/F) wants third degree Taylor
    :returns: f, g
    """
    r_mag = np.linalg.norm(r_vec)

    u = 1 / (r_mag**3)
    z = r_vec @ r_vec_dot / (r_mag**2)
    q = r_vec_dot @ r_vec_dot / (r_mag**2) - u

    f = 1 - tau ** 2 / (2 * r_mag**3) + (r_vec @ r_vec_dot) * tau**3 / (2 * r_mag**5)
    g = tau - tau**3 / (6 * r_mag**3)

    if not want_third:      # if want fourth degree
        f += (3 * u * q - 15 * u * z**2 + u**2) / 24 * tau**4
        g += u * z * tau**4 / 4
    return f, g


def getdEInit(a, e, n, r_vec, r_vec_dot, tau):
    """
    Get initial dE estimate for Newton-Raphson method to determine dE
    """
    if e <= 0.1:
        return n * tau
    r_mag = np.linalg.norm(r_vec)
    sign = r_vec @ r_vec_dot / (n * a**2) * math.cos(n * tau - r_vec @ r_vec_dot / (n * a**2)) + (1 - r_mag / a) \
           * math.sin(n * tau - r_vec @ r_vec_dot / (n * a**2))
    return n * tau + np.sign(sign) * (0.85 * e - r_vec @ r_vec_dot / (n * a**2))


def getdECorr(dE, a, r_vec, r_vec_dot, n, tau):
    """
    Get dE correction for Newton-Raphson method to determine dE
    """
    r_mag = np.linalg.norm(r_vec)
    norm = dE - (1 - r_mag / a) * math.sin(dE) + r_vec @ r_vec_dot / (n * a**2) * (1 - math.cos(dE)) - n * tau
    prime = 1 - (1 - r_mag / a) * math.cos(dE) + r_vec @ r_vec_dot / (n * a**2) * math.sin(dE)
    return - norm / prime


def getdENR(dE, r_vec, r_vec_dot, a, n, tau, threshold):
    """
    Determine dE using Newton-Raphson method
    """
    correction = getdECorr(dE, a, r_vec, r_vec_dot, n, tau)
    counter = 0
    while abs(correction) > threshold:
        dE += correction
        correction = getdECorr(dE, a, r_vec, r_vec_dot, n, tau)
        counter += 1
    return dE


def getFGFunc(tau, r_vec, r_vec_dot, a, e, threshold=1e-12):
    """
    Implements eqs. 109 and 110 functions to get f and g
    """
    r_mag = np.linalg.norm(r_vec)
    n = a**(-3/2)
    dE_i = getdEInit(a, e, n, r_vec, r_vec_dot, tau)
    dE = getdENR(dE_i, r_vec, r_vec_dot, a, n, tau, threshold)
    f = 1 - a / r_mag * (1 - math.cos(dE))
    g = tau + (math.sin(dE) - dE) / n
    return f, g


def getF1G1F3G3(tau1, tau3, r2_vec, r2_dot_vec, flag):
    """
    :param tau1: time at obs 1
    :param tau3: time at obs 3
    :param r2_vec: Earth -> asteroid pos vec
    :param r2_dot_vec: Earth -> asteroid velocity vec
    :param flag: 0: actual equation, 1 : 3rd degree, 2: 4th degree
    :return: f1, g1, f3, g3
    """
    # figure out flag here
    f1 = g1 = f3 = g3 = 0
    if flag == 0:
        a = getA(r2_vec, r2_dot_vec)
        e = getE(getHVec(r2_vec, r2_dot_vec), a)
        f1, g1 = getFGFunc(tau1, r2_vec, r2_dot_vec, a, e)
        f3, g3 = getFGFunc(tau3, r2_vec, r2_dot_vec, a, e)

    elif flag == 1 or flag == 2:
        f1, g1 = getFGTaylor(tau1, r2_vec, r2_dot_vec, flag == 1)
        f3, g3 = getFGTaylor(tau3, r2_vec, r2_dot_vec, flag == 1)
    return f1, g1, f3, g3


def getOrbPlaneR(a, e, E):
    return np.array([a * math.cos(E) - a * e, a * (1 - e**2) ** (1/2) * math.sin(E), 0])


def orbToEcl(orb_el, orb_coords):
    """
    :param orb_el: OrbitalElements object
    :param orb_coords: np array orbital plane position coords
    :return: coordinates transformed to ecliptic plane
    """
    i = orb_el.i
    LOAN = orb_el.LOAN
    APE = orb_el.APE

    peri_trans = np.array([
        [math.cos(APE), -math.sin(APE), 0],
        [math.sin(APE), math.cos(APE), 0],
        [0, 0, 1]
    ])

    incl_trans = np.array([
        [1, 0, 0],
        [0, math.cos(i), - math.sin(i)],
        [0, math.sin(i), math.cos(i)]
    ])

    ecl_x_trans = np.array([
        [math.cos(LOAN), -math.sin(LOAN), 0],
        [math.sin(LOAN), math.cos(LOAN), 0],
        [0, 0, 1]
    ])

    return ecl_x_trans @ incl_trans @ peri_trans @ orb_coords


def eclToEq(ecl_coords, obl=k_obliquity):
    """
    Transform coordinates from ecliptic plane to equatorial plane
    :param ecl_coords: ecliptic coords (np array)
    :param obl: obliquity angle (default = Earth's)
    :return: coordinates transformed to equatorial plane
    """
    return np.array([
        [1, 0, 0],
        [0, math.cos(obl), -math.sin(obl)],
        [0, math.sin(obl), math.cos(obl)]
    ]) @ ecl_coords


def eqToEcl(eq_coords, obl=-k_obliquity):
    """
    Transform coordinates from ecliptic plane to equatorial plane
    :param eq_coords: equatorial coords (np array)
    :param obl: obliquity angle (default = Earth's)
    :return: coordinates transformed to ecliptic plane
    """
    return np.linalg.inv(np.array([
        [1, 0, 0],
        [0, math.cos(obl), -math.sin(obl)],
        [0, math.sin(obl), math.cos(obl)]
    ])) @ eq_coords


def getEarthAst(earth_sun_vec, sun_ast):
    """
    :param earth_sun_vec: E -> S vector
    :param sun_ast: S -> ast vector
    :return: unit vector E -> asteroid
    """
    ret = earth_sun_vec + sun_ast
    norm = np.linalg.norm(ret)
    return ret/norm


def getRADec(earth_ast_vec):
    """
    Get RA,DEC from unit vector from Earth to asteroid
    :param earth_ast_vec: earth asteroid vector
    :return: RA, DEC
    """
    dec = math.asin(earth_ast_vec[2])
    cos_asc = earth_ast_vec[0] / math.cos(dec)
    sin_asc = earth_ast_vec[1] / math.cos(dec)
    asc = getCorrectQuadrant(sin_asc, cos_asc)
    return asc, dec


def getEErr(e, E, M):
    return E - e * math.sin(E) - M


def getEErrPrime(e, E, M):
    return -1 + e * math.cos(E)


def getECorrection(e, E, M):
    return - getEErr(e, E, M) / getEErrPrime(e, E, M)


def getENR(Ei, M, e, k_tol=1e-12, log=False):
    """
    Get Eccentric Anomaly using Newton-Raphson method
    :param Ei: initial guess for E
    :param M: Mean Anomaly
    :param e: eccentricity
    :param k_tol: error tolerance
    :param log: log results
    :return: optimized E
    """
    counter = 0
    error = getECorrection(e, Ei, M)
    while abs(error) > k_tol:
        if log:
            print("Iteration:", counter, "\t\t\t\tCurr E:", Ei)
        Ei -= error

        error = getECorrection(e, Ei, M)
        counter += 1
    if log:
        print("\nTotal Iterations:", counter)
        print("E =", Ei)
        print("Convergence Param:", k_tol)
        print()
    return Ei


def orbElToRADec(orb_el, earth_sun_vec):
    """
    Get RA and DEC from orbital elements
    :param orb_el: OrbitalElements object
    :param earth_sun_vec: pos vector Earth->Sun
    :returns: RA (degrees), DEC (degrees)
    """
    MA = orb_el.MA
    a = orb_el.a
    e = orb_el.e
    E = getENR(MA, MA, e)

    r_vec_orb = getOrbPlaneR(a, e, E)

    r_vec_ecl = orbToEcl(orb_el, r_vec_orb)
    r_vec_eq = eclToEq(r_vec_ecl)

    eq_earth_sun_vec = eclToEq(earth_sun_vec)
    rho = getEarthAst(eq_earth_sun_vec, r_vec_eq)
    ra, dec = getRADec(rho)
    return math.degrees(ra), math.degrees(dec)


def getCorrectObsTime(time_i, dist):
    """
    Correct observation time for lightspeed latency
    :param time_i: initial time
    :param dist: distance traveled (AU)
    :return: corrected time
    """
    return time_i - dist / k_c_AU


def getRhoHatFromObs(ra, dec):
    return np.array([math.cos(ra) * math.sin(dec),
                       math.sin(ra) * math.cos(dec),
                       math.sin(dec)])


def getAllDs(rhos, sun_vecs):
    """
    Use eq. 101b to determine D0, D1J, etc. for 3 observations
    :param rhos: Earth -> asteroid unit vecs (np array)
    :param sun_vecs: Earth -> Sun vecs (np array)
    :return: D vector in form (D0, (D11, D21, D31), (D21, D22...etc.)
    """
    ret = []
    for i in range(3):
        out = np.array(getRowDs(rhos, sun_vecs[i]))
        if len(ret) == 0:
            ret.append(out[0])
        ret.append(out[1:])
    return np.array(ret)


def getRowDs(rhos, sun_vec):
    """
    Use eq. 101b to determine D0, D1J, etc. for observation J
    :param rhos: Earth->asteroid vecs (np array)
    :param sun_vec: Earth -> Sun vec (np array)
    :returns: (D0, D1J, D2J, D3J)
    """
    D_0 = rhos[0] @ (np.cross(rhos[1], rhos[2]))
    D_1J = np.cross(sun_vec, rhos[1]) @ rhos[2]
    D_2J = np.cross(rhos[0], sun_vec) @ rhos[2]
    D_3J = rhos[0] @ np.cross(rhos[1], sun_vec)
    return D_0, D_1J, D_2J, D_3J


def handleMultR2Rho2(r2s, rho2s):
    """
    Query user about which set of (r2, rho2) to use
    """
    if len(r2s) <= 1:
        return r2s[0], rho2s[0]

    print("You have multiple roots!")
    for i in range(len(r2s)):
        print(str(i+1) + ") R2:", r2s[i], ", Rho2", rho2s[i])
    des = input("Please select the one you'd like to use: ")
    while True:
        try:
            des = int(des) - 1
            if des != 1 and des != 2 and des != 3:
                int("truck")        # intentionally go to 'except' circumstance
            break
        except ValueError:
            print("Please enter a valid response!")
    return r2s[des], rho2s[des]


def getCs(f1, f3, g1, g3):
    """
    Implement eqs. 95 and 96 to get c1 and c3
    """
    c1 = g3 / (f1 * g3 - g1 * f3)
    c3 = -g1 / (f1 * g3 - g1 * f3)
    return c1, c3


def getAllRhos(cs, ds):
    """
    Use Scalar Eqs of Range (eq. 101)
    :param cs: (c1, c2, c3)
    :param ds: (d0), (d11, d21, d31), (d12, etc.)
    :return: (rho1, rho2, rho3)
    """
    notD0s = ds[1:]
    rhos = []
    for row in range(3):
        d_row = np.array([notD0s[0][row], notD0s[1][row], notD0s[2][row]])
        rho = d_row @ cs
        rho /= cs[row] * ds[0]
        rhos.append(rho)
    return rhos


def getRs(rhohats, rhomags, sun_vecs):
    """
    Implmeent eq. 2 to get r vectors
    :param rhohats: rhohat vecs (np arrays)
    :param rhomags: rho magnitudes
    :param sun_vecs: Earth -> Sun vecs
    :return: (r1, r2, r3)
    """
    rs = []
    for i in range(min(len(rhohats), len(rhomags), len(sun_vecs))):
        rs.append(rhomags[i] * rhohats[i] - sun_vecs[i])
    return rs


def getR2Dot(f1, f3, g1, g3, rvecs):
    """
    Implements eq. 105 to get initial r2 dot velocity vector
    """
    d1 = -f3 / (f1 * g3 - f3 * g1)
    d3 = f1 / (f1 * g3 - f3 * g1)
    return d1 * rvecs[0] + d3 * rvecs[2]


def getTaus(times):
    tau = k_k * (times[2] - times[0])
    tau1 = k_k * (times[0] - times[1])
    tau3 = k_k * (times[2] - times[1])
    return tau, tau1, tau3


def MOG(obs, want_ecl=False, k_threshold=1e-12, log=False):
    """
    :param obs: Data from 3 observations including:
        - UT Date/Time in seconds from 00:00:00.0 (JD)
        - RA, DEC (sex. hrs., sex. deg.)
        - Earth to Sun vector (AU, eq. cartesian coords., J2000, apparent)
    :param want_ecl: (T/F) wants output to be in ecliptic coordinates
    :param k_threshold: threshold for stopping iterations
    :param log: (T/F) whether or not to log debugging statments
    :returns:
        - Pos/Vel vec for middle obs (ecliptic rect. coords., AU, AU/day)
        - Orbital elements (MA for central observation)
    """
    # SETUP:
    obs_times = [obs[i][0] for i in range(3)]   # JD
    earth_sun_vecs = [obs[i][3] for i in range(3)]
    # earth_sun_vecs = [eclToEq(obs[i][3]) for i in range(3)]

    tau, tau1, tau3 = getTaus(obs_times)
    taus = (tau, tau1, tau3)

    c2 = -1

    rhohats = []
    for observation in obs:
        ra = math.radians(observation[1])
        dec = math.radians(observation[2])
        rhohat = np.array([math.cos(ra) * math.cos(dec),
                           math.sin(ra) * math.cos(dec),
                           math.sin(dec)])
        rhohats.append(rhohat)
    ds = getAllDs(rhohats, earth_sun_vecs)

    # get initial guess for r2 mag
    d2s = np.array([ds[0], ds[1][1], ds[2][1], ds[3][1]])
    r2s_i, rho2s_i = getInitR2sRho2s(taus, earth_sun_vecs[1], rhohats[1], d2s)
    r2_mag, rho2_mag = handleMultR2Rho2(r2s_i, rho2s_i)

    f1, g1 = getInitialFG(tau1, r2_mag)
    f3, g3 = getInitialFG(tau3, r2_mag)

    c1, c3 = getCs(f1, f3, g1, g3)

    rho_mags = getAllRhos(np.array([c1, c2, c3]), ds)
    r_vecs = getRs(rhohats, rho_mags, earth_sun_vecs)

    r2dot_vec = getR2Dot(f1, f3, g1, g3, r_vecs)
    adj_times = [getCorrectObsTime(obs_times[i], rho_mags[i]) for i in range(3)]

    # iteration (do while)
    iters = 0
    fg_flag = 0
    while True:
        tau, tau1, tau3 = getTaus(adj_times)

        f1, g1, f3, g3 = getF1G1F3G3(tau1, tau3, r_vecs[1], r2dot_vec, fg_flag)

        c1, c3 = getCs(f1, f3, g1, g3)
        rho_mags_new = getAllRhos(np.array([c1, -1, c3]), ds)
        r_vecs = getRs(rhohats, rho_mags, earth_sun_vecs)
        r2dot_vec = getR2Dot(f1, f3, g1, g3, r_vecs)
        new_times = [getCorrectObsTime(obs_times[i], rho_mags_new[i]) for i in range(3)]

        diff = np.linalg.norm(rho_mags_new[1] - rho_mags[1])
        if abs(diff) < k_threshold:
            break
        rho_mags = rho_mags_new

        if iters > 200:
            if fg_flag == 2:
                return None
            fg_flag += 1
            iters = 0

        if log:
            print("time correction:", (rho_mags[1] / k_c_AU) * 24 * 3600)
            print("Mag of rho2 change:", diff)
            print("Rho2 mag", rho_mags[1])
            print("Flag:", fg_flag)
            print("Iter:", iters)
        iters += 1

        adj_times = new_times
    if want_ecl:
        ecl_r2_vec = eqToEcl(r_vecs[1], k_obliquity)
        ecl_r2_dot_vec = eqToEcl(r2dot_vec, k_obliquity)
        orb_el = getInfinityStones(adj_times[1], ecl_r2_vec, ecl_r2_dot_vec)
        return ecl_r2_vec, ecl_r2_dot_vec, orb_el
    orb_el = getInfinityStones(adj_times[1], r_vecs[1], r2dot_vec)
    return r_vecs[1], r2dot_vec, orb_el
