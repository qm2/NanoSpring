import random

base = ["A", "T", "C", "G"]


def createRandomGenome(fileName="test", length=5000000):
    """
    Creates a random genome of size length and writes it in fileName.genome
    The first line is an integer representing the length of the genome and the
    second line is the entire genome
    """
    fileName = fileName + ".genome"
    with open(fileName, "w") as f:
        f.write("%d\n" % length)
        f.write(''.join(random.choices(base, k=length)))
        f.write('\n')


def generateRead(genome, readLen, pIn, pDel, pS):
    """
    Generates a read of length readLen from the given genome
    :param genome: The genome from which to generate the reads
    :param pIn: Insertion error rate.
    :param pDel: Deletion error rate.
    :param pS: Substitution error rate.
    :return: (editString, read)
    editString is the string of edits symbolizing how to get from the original
    string to the read string.
    """
    genomeLength = len(genome)
    pos = random.randrange(genomeLength)
    read = ''
    editString = f"{pos}:"
    l = 0
    isUnmodified = False
    nUnmodified = 0
    while (l < readLen):
        r = random.random()

        def writeUnmodified():
            nonlocal isUnmodified
            nonlocal nUnmodified
            nonlocal editString
            if (isUnmodified):
                editString += str(nUnmodified)
                isUnmodified = False
                nUnmodified = 0
        # With probability pIn there is a random insertion
        if (r < pIn):
            baseToAdd = random.choice(base)
            read += baseToAdd
            l += 1
            writeUnmodified()
            editString += "i" + baseToAdd
        # with probability pDel there is a random deletion
        elif (r < pIn + pDel):
            pos += 1
            writeUnmodified()
            editString += "d" + genome[(pos - 1) % genomeLength]
        # With probability pS there is a substitution
        elif (r < pIn + pDel + pS):
            originalBase = genome[pos % genomeLength]
            newBases = base[:]
            newBases.remove(originalBase)
            newBase = random.choice(newBases)
            read += newBase
            l += 1
            pos += 1
            writeUnmodified()
            editString += "s" + originalBase + newBase
        # Otherwise there is just normal transmission
        else:
            read += genome[pos % genomeLength]
            l += 1
            pos += 1
            if not isUnmodified:
                editString += "u"
            isUnmodified = True
            nUnmodified += 1
    writeUnmodified()
    return (editString, read)


def generateReads(fileName="test", coverage=25, readLen=10000, pIn=0.00, pDel=0.00, pS=0.05):
    """
    Generate random reads of length readLen with given coverage from the genome
    in the file fileName.genome and stores the reads in fileName.reads. Each
    read is stored in two lines:
    readPos:edit string
    readData
    readPos is the position in the original genome (index starting from 0) where
    the read starts
    edit string is the string of edits symbolizing how to get from the original
    string to the read string.
    Example:
    u20iAdCsTG
    20 Unmodified Insert A, Delete C, Substitute G for T (T->G).
    Each read has fixed length readLen and if the end of the genome is reached,
    the read sequence wraps around to the beginning.
    :param pIn: Insertion error rate.
    :param pDel: Deletion error rate.
    :param pS: Substitution error rate.
    """
    genomeFileName = fileName + ".genome"
    genome = ''
    genomeLength = 0
    with open(genomeFileName, "r") as f:
        genomeLength = int(f.readline())
        genome = f.readline()[:-1]
    readFileName = fileName + ".reads"
    with open(readFileName, "w") as f:
        numReads = int(genomeLength * coverage / readLen)
        for i in range(numReads):
            (editString, read) = generateRead(genome, readLen, pIn, pDel, pS)
            f.write(editString)
            f.write('\n')
            f.write(read)
            f.write('\n')


if __name__ == '__main__':
    createRandomGenome(length=100000)
    generateReads(readLen=10000)
