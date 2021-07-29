class OrbitalElements:
    def __init__(self, a, e, i, LOAN, APE, MA, PET):
        self.a = a
        self.e = e
        self.i = i
        self.LOAN = LOAN
        self.APE = APE
        self.MA = MA
        self.PET = PET


    def __str__(self):
        return "A: " + str(self.a) + ", E: " + str(self.e) + ", I (rad): " + str(self.i)\
               + ", Long. of Asc. Node (rad): " + str(self.LOAN) + ", Arg. of Perihelion (rad): " + str(self.APE)\
               + ", Mean Anomaly (rad): " + str(self.MA) + ", Time of Peri. Passage (JD): " + str(self.PET)
