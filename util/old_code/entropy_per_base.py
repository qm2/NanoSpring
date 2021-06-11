import numpy as np


def entropy_per_base(ps, pi, pd):
    def xlogx(p):
        if p == 0:
            return 0
        return p*np.log2(p)
    return (-xlogx(1-ps-pd-pi) - 3*xlogx(ps/3) - xlogx(pd) - 4*xlogx(pi/4))/(1-pd)


if __name__ == '__main__':
    import sys
    if (len(sys.argv) < 4):
        print("Usage ./entropy_per_base.py ps pi pd")
    else:
        print(entropy_per_base(
            float(sys.argv[1]),
            float(sys.argv[2]),
            float(sys.argv[3]),
        ))
