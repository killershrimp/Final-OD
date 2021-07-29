# AUTHOR: Alex Yun
# ASSIGNMENT: Final OD
# COLLABORATORS: Justin J., Emma C., Neil S., Saraliba A., Daisy Z.

# STATUS: runs fine, scalars same as test output to ~5 decimal places


from libs.odlib import *


obs = parseODData("sunnyinp.txt")
mog_ret = MOG(obs, log=True, want_ecl=True)
if mog_ret is None:
    print("Didn't converge!")
    exit(1)
pos_vec = mog_ret[0]
vel_vec = mog_ret[1]
orb_el_JD = mog_ret[2]

des_JD = getJD(2021, 7, 8) + timeToJD("00:22:04")
orb_el_JD = incrementMA(orb_el_JD, des_JD)

print("Ecliptic Pos Vec:", pos_vec)
print("Ecliptic Vel Vec:", vel_vec)
print("Range to Asteroid for Middle Observation (elc):", np.linalg.norm(pos_vec), "AU")
print("Orbital Elements:")
print("\ta (AU):", orb_el_JD[0].a)
print("\te:", orb_el_JD[0].e)
print("\ti (deg):", math.degrees(orb_el_JD[0].i))
print("\tArg. of Perihelion (deg):", math.degrees(orb_el_JD[0].APE))
print("\tLong. of Asc. Node (deg):", math.degrees(orb_el_JD[0].LOAN))
print("\tMean Anomaly for 07/08/2021 at 00:22:04 UT (deg):", math.degrees(orb_el_JD[0].MA))

ra, dec = orbElToRADec(orb_el_JD[0], obs[1][3])
print("Expected RA at 07/08/2021 at 00:22:04 UT:", decDegToHr(ra))
print("Expected Dec at 07/08/2021 at 00:22:04 UT:", decDegToSexDeg(dec))

