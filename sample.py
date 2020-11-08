import random


def sample(oldfile, newfile, p):
    with open(oldfile, "r") as oldf:
        with open(newfile, "w") as newf:
            while True:
                line = oldf.readline()
                if not line:
                    break
                r = random.random()
                if r < p:
                    for i in range(3):
                        newf.writelines([line])
                        line = oldf.readline()
                    newf.writelines([line])
                else:
                    for i in range(3):
                        line = oldf.readline()


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 4:
        print("Usage: ./sample.py old new p")
    else:
        sample(
            sys.argv[1],
            sys.argv[2],
            float(sys.argv[3]),
        )
