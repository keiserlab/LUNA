
# class OpenEyeChargeModel():

    # def __init__():

        # self.

def perceive_formal_charge(rdAtom):
    currentFormalCharge = rdAtom.GetFormalCharge()
    formalCharge = None

    if (currentFormalCharge == 0):
        atomicNum = rdAtom.GetAtomicNum()
        valence = rdAtom.GetExplicitValence()

        # Hydrogen
        if atomicNum == 1:
            if valence != 1:
                formalCharge = 1
        # Carbon
        elif atomicNum == 6:
            if valence == 3:
                hasPolarNeighbor = False
                for neighbor in rdAtom.GetNeighbors():
                    # Polar neighbor: N, O or S
                    if neighbor.GetAtomicNum() in (7, 8, 16):
                        hasPolarNeighbor = True
                        break
                if hasPolarNeighbor is True:
                    formalCharge = 1
                else:
                    formalCharge = -1
        # Nitrogen
        elif atomicNum == 7:
            if valence == 2:
                formalCharge = -1
            elif valence == 4:
                formalCharge = 1
        # Oxygen
        elif atomicNum == 8:
            if valence == 1:
                formalCharge = -1
            elif valence == 3:
                formalCharge = 1
        # Phosphorus
        elif atomicNum == 15:
            if valence == 4:
                formalCharge = 1
        # Sulfur
        elif atomicNum == 16:
            if valence == 1:
                formalCharge = -1
            elif valence == 3:
                formalCharge = 1
            elif valence == 5:
                formalCharge = -1
            elif valence == 4:
                degree = rdAtom.GetDegree()
                if (degree == 4):
                    formalCharge = 2
        # Chlorine
        elif atomicNum == 17:
            if valence == 0:
                formalCharge = -1
            if valence == 4:
                formalCharge = 3
        # Fluorine, Bromine, Iodine
        elif atomicNum == 9 or atomicNum == 35 or atomicNum == 53:
            if valence == 0:
                formalCharge = -1
        # Magnesium, Calcium, Zinc
        elif atomicNum == 12 or atomicNum == 20 or atomicNum == 30:
            if valence == 0:
                formalCharge = 2
        # Lithium, Sodium, Potassium
        elif atomicNum == 3 or atomicNum == 11 or atomicNum == 19:
            if valence == 0:
                formalCharge = 1
        # Boron
        # If the valence is four, the formal charge is +1.
        elif atomicNum == 5:
            if valence == 4:
                formalCharge = 1
        else:
            # Partial charge > 0.4 => +1
            # Partial charge < 0.4 => -1
            # Consider bins or rounding


    return formalCharge
    # Iron
    # If the valence is zero, the formal charge is +3 if the partial charge is 3.0, and +2 otherwise.
    # elif atomicNum == 26:
        # if valence == 0:
            # formalCharge = 3
    # Copper
    # If the valence is zero, the formal charge is +2 if the partial charge is 2.0, and +1 otherwise.





            # Sulfur
                # If the valence is 1, the formal charge is -1,
                # if the valence is three the formal charge is +1,
                # if the valence is 5, the formal charge is -1,
                # if the valence is four and the degree is four the charge is +2.


            # Chlorine
                # If the valence is 0 the formal charge is -1, if the valence is four the formal charge is +3.

            # Fluorine, Bromine, Iodine
            # If the valence is zero, the formal charge is -1.

            # Magnesium, Calcium, Zinc
            # If the valence is zero, the formal charge is +2.
            # Lithium, Sodium, Potassium
            # If the valence is zero, the formal charge is +1.
# Iron
# If the valence is zero, the formal charge is +3 if the partial charge is 3.0, and +2 otherwise.
# Copper
# If the valence is zero, the formal charge is +2 if the partial charge is 2.0, and +1 otherwise.
