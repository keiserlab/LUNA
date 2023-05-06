import argparse


CHARGE_METHODS = ['am1-bcc', 'gasteiger']


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', dest="input", type=str, required=True,
                        help="an input molecular file.")

    parser.add_argument('-o', dest="output", type=str, required=True,
                        help="the output file.")

    parser.add_argument('-nc', dest="total_charge", type=int, required=True,
                        help="net charge.")

    parser.add_argument('-c', dest="charge_method", type=str,
                        default="am1-bcc",
                        help="Charge method is either 'am1-bcc' or "
                             "'gasteiger'. The default value is 'am1-bcc'.")

    args = parser.parse_args()

    if args.charge_method not in CHARGE_METHODS:
        raise ValueError("Unknown charge method: '%s'." % args.charge_method)

    from chimera import runCommand as rc
    from AddCharge import addNonstandardResCharges
    from chimera.selection import currentResidues

    rc("open %s" % args.input)

    rc("sel #0")  # select the ligand
    ligand = currentResidues()[0]

    addNonstandardResCharges([ligand], args.total_charge,
                             method=args.charge_method)

    rc("write format mol2 #0 %s" % args.output)

    chimera.closeSession()

    rc('stop now')


if __name__ == "__main__":
    main()
