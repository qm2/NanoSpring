import sys
bases = ['A', 'T', 'C', 'G']

def toBaseString(n):
    result = []
    while (n > 0):
        result.insert(0, bases[n % 4])
        n = n // 4
    return(''.join(result))

if __name__ == "__main__":
    # execute only if run as a script
    print(toBaseString(int(sys.argv[1])))
