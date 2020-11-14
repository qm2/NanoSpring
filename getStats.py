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


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: ./getStatst.py file")
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
