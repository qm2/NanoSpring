def getStats(f):
    types = 'MDISH'
    count = {}
    for t in types:
        count[t] = 0
    while(True):
        line = f.readline()
        if not line:
            break
        if len(line) == 0 or line[0] == '@':
            continue
        fields = line.split('\t')
        cigar = fields[5]
        import re
        edits = re.split('([{}])'.format(types), cigar)[:-1]
        num = len(edits)
        for i in range(0, num, 2):
            count[edits[i+1]] += int(edits[i])
    return count


def cigar_main():
    import sys
    if len(sys.argv) < 2:
        print("Usage: ./getStats.py file")
    else:
        with open(sys.argv[1]) as f:
            count = getStats(f)
            print(count)
            numBasesRead = 0
            # for (_, c) in count.items():
            # numBasesRead += c
            for t in 'MID':
                numBasesRead += count[t]
            numIns = count['I']
            numDel = count['D']
            print(
                f"numBasesRead:{numBasesRead} numIns:{numIns} numDel:{numDel}")
            print(
                f"pIns:{numIns / numBasesRead} pDel:{numDel / numBasesRead}"
            )


def getCSStats(f):
    count = {'u': 0, 'i': 0, 'd': 0, 's': 0}
    while(True):
        line = f.readline()
        if not line:
            break
        if len(line) == 0 or line[0] == '@':
            continue
        fields = line.split('\t')
        cs = ''
        for field in fields:
            if len(field) < 5:
                continue
            if field[0:5] == "cs:Z:":
                cs = field[5:]
        import re
        edits = re.split('([:+*-])', cs)[1:]
        num = len(edits)
        for i in range(0, num, 2):
            if edits[i] == ':':
                count['u'] += int(edits[i+1])
            elif edits[i] == '-':
                count['d'] += len(edits[i+1])
            elif edits[i] == '+':
                count['i'] += len(edits[i+1])
            elif edits[i] == '*':
                count['s'] += 1
    return count


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: ./getStats.py file")
    else:
        with open(sys.argv[1]) as f:
            count = getCSStats(f)
            print(count)
            numBasesRead = count['u'] + count['i'] + count['s'] + count['d']
            print(f"pS:{count['s'] / numBasesRead}")
            print(f"pIns:{count['i'] / numBasesRead}")
            print(f"pDel:{count['d'] / numBasesRead}")
