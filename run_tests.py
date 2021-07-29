from tests.odlib_tests import *


k_want_logs = False

# ODLIB ANGLE CONVERSION TESTS
runTestHrsToDecDeg(log=k_want_logs)
runTestSexDegToDecDeg(log=k_want_logs)
runTestDecDegToHr(log=k_want_logs)
runTestDecDegToSexDeg(log=k_want_logs)

# OD LIB TESTS:
runTestGetInfinityStones(log=k_want_logs)
runTestGetInitR2Rho2(log=k_want_logs)
runTestGetF1G1F3G3(log=k_want_logs)
runTestOrbElToRADec(log=k_want_logs)
