"""
Mix of unit tests and integration tests. Tests run in "run_tests.py"
"""

from libs.odlib import *
import numpy as np


# UNIT CONVERSION TESTS
def testSexDegToDecDeg(inp, des_out, k_threshold=1e-4, log=False):
    der_out = sexDegToDecDeg(inp)
    works = inRange(der_out, des_out, k_threshold)
    if log:
        print("\t", "PASSED" if works else "FAILED", "test for DEC DEG")
        print("\t\tDesired Value:", des_out)
        print("\t\tDerived Value:", der_out)
    return works


def runTestSexDegToDecDeg(log=False):
    print("TESTING sexDegToDecDeg")
    inp = ["18:18:18.18", "0:0:0"]
    des_out = [18.30505000, 0]

    all_works = True
    for i in range(len(inp)):
        works = testSexDegToDecDeg(inp[i], des_out[i], log=log)
        print("\tTest", i+1, "PASSED" if works else "FAILED")
        if not all_works or not works:
            all_works = False

    if all_works:
        print("ALL TESTS FOR sexDegToDecDeg PASSED")
    else:
        print("SOME TESTS FOR sexDegToDecDeg FAILED")
    print()


def testHrsToDecDeg(inp, des_out, k_threshold=1e-6, log=False):
    der_out = hrsToDecDeg(inp)
    works = inRange(der_out, des_out, k_threshold)
    if log:
        print("\t", "PASSED" if works else "FAILED", "test for DEC DEG")
        print("\t\tDesired Value:", des_out)
        print("\t\tDerived Value:", der_out)
    return works


def runTestHrsToDecDeg(log=False):
    print("TESTING hrsToDecDeg")
    inp = ["18:18:18.18", "0:0:0"]
    des_out = [274.57575000, 0]

    all_works = True
    for i in range(len(inp)):
        works = testHrsToDecDeg(inp[i], des_out[i], log=log)
        print("\tTest", i+1, "PASSED" if works else "FAILED")
        if not all_works or not works:
            all_works = False

    if all_works:
        print("ALL TESTS FOR hrsToDecDeg PASSED")
    else:
        print("SOME TESTS FOR hrsToDecDeg FAILED")
    print()


def testDecDegToHr(inp, des_out, log=False):
    des_angle = inp[0]
    round_places = inp[1]
    der_out = decDegToHr(des_angle, round_places=round_places)
    works = der_out == des_out
    if log:
        print("\t", "PASSED" if works else "FAILED", "test for SEX DEC conversion")
        print("\t\tDesired Value:", des_out)
        print("\t\tDerived Value:", der_out)
    return works


def runTestDecDegToHr(log=False):
    print("TESTING decDegToHr")
    inp = [(0, 0), (274.5750, 0)]
    des_out = ["00:00:00.0", "18:18:18.0"]

    all_works = True
    for i in range(len(inp)):
        works = testDecDegToHr(inp[i], des_out[i], log=log)
        print("\tTest", i+1, "PASSED" if works else "FAILED")
        if not all_works or not works:
            all_works = False

    if all_works:
        print("ALL TESTS FOR decDegToHr PASSED")
    else:
        print("SOME TESTS FOR decDegToHr FAILED")
    print()


def testDecDegToSexDeg(inp, des_out, log=False):
    des_angle = inp[0]
    round_places = inp[1]
    der_out = decDegToSexDeg(des_angle, round_places=round_places)
    works = der_out == des_out
    if log:
        print("\t", "PASSED" if works else "FAILED", "test for HR conversion")
        print("\t\tDesired Value:", des_out)
        print("\t\tDerived Value:", der_out)
    return works


def runTestDecDegToSexDeg(log=False):
    print("TESTING decDegToSexDeg")
    inp = [(0, 0), (18.30500000, 0), (-18.30500000, 0)]
    des_out = ["+00:00:00.0", "+18:18:18.0", "-18:18:18.0"]

    all_works = True
    for i in range(len(inp)):
        works = testDecDegToSexDeg(inp[i], des_out[i], log=log)
        print("\tTest", i+1, "PASSED" if works else "FAILED")
        if not all_works or not works:
            all_works = False

    if all_works:
        print("ALL TESTS FOR decDegToSexDeg PASSED")
    else:
        print("SOME TESTS FOR decDegToSexDeg FAILED")


# ORBITAL DETERMINATION TESTS
def testGetH(inp, des_out, print_out=True):
    """
    Tests
    :param inp: input for test case
    :param des_out: desired output for test case
    :return: whether test case passed, derived output
    """
    k_epsilon = 1e-5
    der_out = getHVec(inp[0], inp[1])
    test_pass = inRange(des_out[0], der_out[0], k_epsilon) \
                and inRange(des_out[1], der_out[1], k_epsilon) \
                and inRange(des_out[2], der_out[2], k_epsilon)

    if print_out:
        print("PASSED" if test_pass else "FAILED", "ANGULAR MOMENTUM TEST for:")
        print("\tINPUT", inp)
        print("\tDESIRED OUTPUT", des_out)
        if not test_pass:
            print("\tDerived Output", der_out)
    return test_pass


def testGetInfinityStones(JD, inp_state, des_stones, log=False):
    """
    :param JD: test case JD
    :param inp_state: (pos_vec, vel_vec)
    :param des_stones: tuple of desired orbital elements
    :param log: whether or not want to print errors
    :return: T/F test passed
    """
    # semi-major axis, eccentricity, inclination, long. of asc. node, arg. of peri, mean anomaly, last peri passage time
    stone_order = ["a", "e", "i", "LOAN", "APE", "MA", "PET"]
    der_stones = getInfinityStones(JD, inp_state[0], inp_state[1], want_tup=True)[0]  # don't test JD
    k_acceptable_percent = 0.02 / 100
    works = True
    for element in range(len(der_stones)):
        des = des_stones[element]
        der = der_stones[element]
        this_works = inRange(des, der, k_acceptable_percent * des)
        if log:
            print("\t\tTEST", "PASSED" if this_works else "FAILED", "for element", stone_order[element], ":")
            print("\t\t\tDesired:", des)
            print("\t\t\tDerived:", der)
        if not this_works or not works:
            works = False
    print("\tPASSED" if works else "FAILED", "TEST FOR getInfinityStones")
    return works


def runTestGetInfinityStones(log=False):
    # JD, pos, vel
    inp_JD_pos_vel = [[2458312.5, [0.38564181, -1.22768327, 0.47096453], [0.6652513, 0.1476058, 0.2211997]],
        [2458312.541666667, [0.3861186, -1.22757735, 0.47112302], [0.66514382, 0.14794774, 0.22106849]],
        [2458312.583333333, [0.38659531, -1.22747118, 0.47128143], [0.66503622, 0.14828959, 0.22093726]],
        [2458312.625, [0.38707193, -1.22736477, 0.47143974], [0.66492851, 0.14863133, 0.22080602]],
        [2458312.666666667, [0.38754849, -1.22725812, 0.47159796], [0.66482069, 0.14897297, 0.22067476]],
        [2458312.708333333, [0.38802496, -1.22715122, 0.47175608], [0.66471276, 0.14931451, 0.22054348]],
        [2458312.75, [0.38850136, -1.22704407, 0.47191411], [0.66460472, 0.14965596, 0.22041218]],
        [2458312.791666667, [0.38897768, -1.22693668, 0.47207204], [0.66449657, 0.1499973, 0.22028087]],
        [2458312.833333333, [0.38945392, -1.22682905, 0.47222988], [0.66438831, 0.15033854, 0.22014954]],
        [2458312.875, [0.38993009, -1.22672117, 0.47238763], [0.66427993, 0.15067969, 0.22001819]],
        [2458312.916666667, [0.39040617, -1.22661305, 0.47254528], [0.66417145, 0.15102073, 0.21988682]],
        [2458312.958333333, [0.39088218, -1.22650468, 0.47270284], [0.66406286, 0.15136168, 0.21975544]],
        [2458313.0, [0.39135812, -1.22639607, 0.4728603], [0.66395416, 0.15170252, 0.21962404]],
        [2458313.041666667, [0.39183397, -1.22628722, 0.47301767], [0.66384535, 0.15204327, 0.21949263]],
        [2458313.083333333, [0.39230974, -1.22617812, 0.47317495], [0.66373642, 0.15238392, 0.21936119]],
        [2458313.125, [0.39278544, -1.22606877, 0.47333213], [0.66362739, 0.15272447, 0.21922974]],
        [2458313.166666667, [0.39326106, -1.22595918, 0.47348922], [0.66351825, 0.15306491, 0.21909827]],
        [2458313.208333333, [0.3937366, -1.22584935, 0.47364621], [0.66340899, 0.15340526, 0.21896679]],
        [2458313.25, [0.39421206, -1.22573928, 0.47380311], [0.66329963, 0.15374551, 0.21883529]],
        [2458313.291666667, [0.39468745, -1.22562896, 0.47395991], [0.66319016, 0.15408566, 0.21870377]],
        [2458313.333333333, [0.39516275, -1.22551839, 0.47411662], [0.66308058, 0.15442572, 0.21857223]]]

    # JD, a, e, math.radians(i), math.radians(LOAN), math.radians(APE), math.radians(MA), PET
    inp_JD_orb_el = [[2458312.5, 1.056800391014609, 0.3442329497516519, 0.43904206003436086, 4.123130888314789, 4.45939770630494,
      2.434945841158005, 2458158.7209008667],
     [2458312.541666667, 1.056800376888049, 0.3442329564821012, 0.43904205603684243, 4.123130875387333,
      4.459397667364686, 2.4356056754799806, 2458158.7208987554],
     [2458312.583333333, 1.05680036277521, 0.344232963205408, 0.4390420520453394, 4.123130862469327, 4.459397628429943,
      2.436265509818305, 2458158.720896642],
     [2458312.625, 1.056800348676084, 0.3442329699215829, 0.43904204805984987, 4.123130849560788, 4.459397589500658,
      2.436925344172978, 2458158.7208945258],
     [2458312.666666667, 1.05680033459066, 0.3442329766306361, 0.43904204408037195, 4.1231308366617325,
      4.4593975505767895, 2.4375851785439995, 2458158.7208924075],
     [2458312.708333333, 1.056800320518931, 0.3442329833325784, 0.43904204010690406, 4.123130823772176,
      4.459397511658289, 2.438245012931376, 2458158.720890287],
     [2458312.75, 1.056800306460885, 0.3442329900274205, 0.4390420361394445, 4.123130810892135, 4.459397472745109,
      2.4389048473351123, 2458158.7208881634],
     [2458312.791666667, 1.056800292416515, 0.3442329967151722, 0.4390420321779911, 4.123130798021624,
      4.459397433837206, 2.4395646817552024, 2458158.720886038],
     [2458312.833333333, 1.056800278385811, 0.344233003395844, 0.43904202822254235, 4.123130785160657, 4.45939739493453,
      2.440224516191658, 2458158.7208839105],
     [2458312.875, 1.056800264368763, 0.3442330100694472, 0.4390420242730961, 4.123130772309251, 4.459397356037038,
      2.440884350644474, 2458158.72088178],
     [2458312.916666667, 1.056800250365362, 0.3442330167359917, 0.4390420203296508, 4.123130759467423,
      4.459397317144683, 2.441544185113653, 2458158.7208796474],
     [2458312.958333333, 1.056800236375599, 0.3442330233954877, 0.43904201639220425, 4.1231307466351845,
      4.459397278257415, 2.4422040195992065, 2458158.720877513],
     [2458313.0, 1.056800222399465, 0.3442330300479461, 0.43904201246075486, 4.123130733812554, 4.4593972393751935,
      2.442863854101131, 2458158.7208753754],
     [2458313.041666667, 1.05680020843695, 0.3442330366933772, 0.4390420085353003, 4.123130720999545, 4.459397200497971,
      2.443523688619425, 2458158.7208732357],
     [2458313.083333333, 1.056800194488045, 0.3442330433317909, 0.43904200461583875, 4.1231307081961726,
      4.4593971616256995, 2.4441835231541007, 2458158.720871094],
     [2458313.125, 1.05680018055274, 0.3442330499631991, 0.4390420007023684, 4.1231306954024545, 4.459397122758333,
      2.444843357705155, 2458158.7208689493],
     [2458313.166666667, 1.056800166631027, 0.3442330565876111, 0.4390419967948873, 4.123130682618402,
      4.459397083895829, 2.4455031922725854, 2458158.7208668026],
     [2458313.208333333, 1.056800152722896, 0.3442330632050378, 0.43904199289339313, 4.123130669844032,
      4.45939704503814, 2.4461630268564076, 2458158.7208646536],
     [2458313.25, 1.056800138828338, 0.3442330698154896, 0.4390419889978842, 4.12313065707936, 4.459397006185217,
      2.446822861456615, 2458158.720862502],
     [2458313.291666667, 1.056800124947343, 0.3442330764189768, 0.4390419851083583, 4.123130644324401,
      4.459396967337019, 2.447482696073212, 2458158.720860348],
     [2458313.333333333, 1.056800111079902, 0.3442330830155106, 0.43904198122481347, 4.123130631579168,
      4.459396928493498, 2.4481425307062032, 2458158.720858192]]

    all_pass = True

    print("TESTING getInfinityStones")
    for i in range(len(inp_JD_pos_vel)):
        sample = inp_JD_pos_vel[i]
        des_out = inp_JD_orb_el[i][1:]  # get rid of JD first element
        JD = sample[0]
        pos_vec = np.array(sample[1])
        vel_vec = np.array(sample[2])

        if log:
            print("\nTest case JD:", JD)
        inp_state = (pos_vec, vel_vec)
        test_pass = testGetInfinityStones(JD, inp_state, des_out, log=log)
        if not (test_pass and all_pass):
            all_pass = False


def testGetInitR2Rho2(inp, des_outp, log=False):
    """
    Tests getInitR2Rho2 method
    :param inp: inputs for test cases
    :param des_outp: desired outputs based on inp
    :param log: T/F want log
    :return: T/F test worked
    """
    der_outp = getInitR2sRho2s(inp[0], inp[1], inp[2], inp[3])
    works = True
    for i in range(min(len(der_outp[0]), len(der_outp[1]))):
        r2_works = inRange(des_outp[1][i], der_outp[1][i], 1e-6)

        if log:
            print("\t\tTEST", "PASSED" if r2_works else "FAILED", "for r2")
            print("\t\t\tDesired r2:", des_outp[0][i])
            print("\t\t\tDerived r2:", der_outp[0][i])

        rho2_works = inRange(des_outp[0][i], der_outp[0][i], 1e-6)
        if log:
            print("\t\tTEST", "PASSED" if rho2_works else "FAILED", "for rho2")
            print("\t\t\tDesired rho2:", des_outp[1][i])
            print("\t\t\tDerived rho2:", der_outp[1][i])
        if not rho2_works or not r2_works or not works:
            works = False
    return works


def runTestGetInitR2Rho2(log=False):
    print("TESTING getInitR2Rho2:")
    taus1 = np.array([-0.15481889055, 0.15481889055, 0.3096377811])
    sun21 = np.array([-0.2398478458274071, 0.9065739917845802, 0.3929623749770952])
    rhohat21 = np.array([-0.8518563498182248, -0.2484702599212149, 0.4610892421311239])
    Ds1 = np.array([-0.0010461861084885213, -0.17297581974209159, -0.17201260125558127, -0.16712421570714076])
    inp1 = [taus1, sun21, rhohat21, Ds1]

    roots1 = np.array([1.014886625023963, 1.2932503440012362, 1.5851855408957922])
    rhos1 = np.array([0.012425430826237482, 0.9753289007273918, 1.386900042701193])
    des_outp1 = [roots1, rhos1]
    print("\tTest 1: PASSED" if testGetInitR2Rho2(inp1, des_outp1, log) else "\tTest 2: FAILED")

    taus2 = np.array([-0.1720209895, 0.1720209895, 0.344041979])
    sun22 = np.array([-0.2683394448727136, 0.8997620807182745, 0.3900022331276332])
    rhohat22 = np.array([0.052719013914983195, -0.9121555187306237, 0.40643943610469035])
    Ds2 = np.array([0.0011797570297704812, 0.052586743761143424, 0.05848153743706686, 0.06274019190783499])
    inp2 = [taus2, sun22, rhohat22, Ds2]

    roots2 = np.array([0.7775650888697365, 1.0071727095060097, 1.266484861890734])
    rhos2 = np.array([-0.8448416611523188, -0.01440831314501545, 0.33742943029525474])
    des_outp2 = [roots2, rhos2]
    print("\tTest 2: PASSED" if testGetInitR2Rho2(inp2, des_outp2, log) else "\tTest 2: FAILED")


def testGetF1G1F3G3(inp, des_outp, flag, log=False):
    k_epsilon = 1.1e-6
    taus = inp[0]
    r2 = inp[1]
    r2_dot = inp[2]

    der_outp = getF1G1F3G3(taus[0], taus[1], r2, r2_dot, flag)

    k_test_order = ["f1", "g1", "f3", "g3"]
    all_work = True
    for i in range(4):
        works = inRange(der_outp[i], des_outp[i], k_epsilon)
        if log:
            print("\t", "PASSED" if works else "FAILED", "test for", k_test_order[i])
            print("\t\tDesired Value:", des_outp[i])
            print("\t\tDerived Value:", der_outp[i])
        all_work = works and all_work
    return all_work


def runTestGetF1G1F3G3(log=False):
    print("TESTING getf1g1f3g3:")

    # TEST 1
    tau1 = -0.32618569435308475
    tau3 = 0.050840808143482484
    r2 = np.array([0.26640998194891174, -1.382856212643199, -0.505199925482389])
    r2dot = np.array([0.8439832722802604, -0.39937767878456487, 0.14200790188593015])
    f1 = 0.9823546285441103
    f3 = 0.9996202204156108
    g1 = -0.3241600025463262
    g3 = 0.050834422945470005
    inp1 = ((tau1, tau3), r2, r2dot)
    des_outp1 = (f1, g1, f3, g3)
    test1 = testGetF1G1F3G3(inp1, des_outp1, -1, log)
    print("\tTest 1", "PASSED" if test1 else "FAILED")

    # TEST 2
    tau1 = -0.32618617484601165
    tau3 = 0.0508408854033231
    r2 = np.array([0.26799552002875776, - 1.3726277901924608, - 0.5026729612047128])
    r2dot = np.array([0.8456809141954584, - 0.3838382184712308, 0.14215854191172816])
    f1 = 0.9821596284506794
    f3 = 0.9996124342607986
    g1 = -0.32442392608030396
    g3 = 0.05083421257607972
    inp2 = ((tau1, tau3), r2, r2dot)
    des_outp2 = (f1, g1, f3, g3)
    test2 = testGetF1G1F3G3(inp2, des_outp2, 3, log)
    print("\tTest 2", "PASSED" if test2 else "FAILED")

    # TEST 3
    tau1 = -0.3261857571141891
    tau3 = 0.05084081855693949
    r2 = np.array([0.26662393644794813, - 1.381475976476564, - 0.5048589337503169])
    r2dot = np.array([0.8442117090940343, - 0.39728396707075087, 0.14202728258915864])
    f1 = 0.9823149707782799
    f3 = 0.99961917185922
    g1 = -0.32418770657924106
    g3 = 0.05083441832100904
    inp3 = ((tau1, tau3), r2, r2dot)
    des_outp3 = (f1, g1, f3, g3)
    test3 = testGetF1G1F3G3(inp3, des_outp3, 4, log)
    print("\tTest 3", "PASSED" if test3 else "FAILED")
    if test1 and test2 and test3:
        print("getf1g1f3g3: ALL TESTS PASSED")
    else:
        print("getf1g1f3g3: SOME TESTS FAILED")


def testOrbElToRADec(inp, des_out, k_threshold=1e-2, log=False):
    test_orb_el = inp[0]
    earth_sun_vec = inp[1]
    der_out = orbElToRADec(test_orb_el, earth_sun_vec)

    k_test_order = ["RA", "DEC"]
    all_work = True
    for i in range(len(k_test_order)):
        works = inRange(der_out[i], des_out[i], k_threshold)
        if log:
            print("\t", "PASSED" if works else "FAILED", "test for", k_test_order[i])
            print("\t\tDesired Value:", des_out[i])
            print("\t\tDerived Value:", der_out[i])
            print("\t\tPercent Error:", abs((der_out[i] - des_out[i]) / des_out[i]))
        all_work = works and all_work
    return all_work


def runTestOrbElToRADec(log=False):
    orb_el = OrbitalElements(
        a=1.056800055745494E+00,
        e=3.442331093323161E-01,
        i=math.radians(2.515525166662531E+01),
        LOAN=math.radians(2.362379806551942E+02),
        APE=math.radians(2.555046142766286E+02),
        MA=math.radians(1.404194576239256E+02),
        PET=2458158.720849543810
    )
    JD = 2458313.500000000
    earth_sun_vec = (2458333.5,
                     np.array([-6.57401119e-01,  7.73019293e-01, -3.33404010e-05]),
                     np.array([-7.45747018e-01, -6.43813191e-01,  2.17135130e-05]))
    JPL_RA_Dec = (hrsToDecDeg("17:42:20.42"), sexDegToDecDeg("31:52:31.3"))

    des_JD = getJD(2018, 8, 3)

    orb_el, JD = incrementMA((orb_el, JD), des_JD)
    print("TESTING orbElToRADec:")
    test1 = testOrbElToRADec((orb_el, earth_sun_vec[1]), JPL_RA_Dec, log=log)
    print("\tTest 1", "PASSED" if test1 else "FAILED")
    if test1:
        print("orbElToRADec: ALL TESTS PASSED")
    else:
        print("orbElToRADec: SOME TESTS FAILED")



def testGetCorrectQuadrant(log=False):
    """
    :return: list of failed getCorrectQuadrant tests
    """
    k_epsilon = 1e-6
    test_angles = [math.pi, 7/8 * math.pi, 20/12 * math.pi, -math.pi / 2, 11/7 * math.pi]
    failed = []

    for i in range(len(test_angles)):
        angle = test_angles[i]
        des_angle = bound0To2Pi(angle)
        der_angle = getCorrectQuadrant(math.sin(angle), math.cos(angle))
        passed = inRange(des_angle, der_angle, k_epsilon)
        if not passed:
            failed.append(des_angle)
        if log:
            print("PASSED" if passed else "FAILED", "test for " + str(i) + ", angle:", angle)
            print("\tDesired angle:", des_angle)
            print("\tDerived angle", der_angle)
    return failed



# TODO add integration test for MOG
