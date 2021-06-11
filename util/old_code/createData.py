import scipy.stats as st
import unittest
import random
import sys
import os
import numpy as np
import configparser


def draw_const_dis(mean, size=1, genomeLen=np.inf):
    return np.full(size, mean)

#------------------------------------------------------------------------------#
# Taken and modified from the code of DeepSimulator


def draw_beta_dis(mean, size=1, genomeLen=np.inf):
    """
    beta distribution with parameters
    1.7780912583853339, 
    7.8927794827841975, 
    316.75817723758286, 
    34191.25716704056
    :param mean: the rough mean of the read length
    :param size: the number of output sequence
    """
    samples = st.beta.rvs(1.778, 7.892, 316.758, 34191.257, size=size)
    samples = samples * mean / 6615.0
    samples = samples.astype(int)
    samples = np.clip(samples, 1, genomeLen)
    return samples


def draw_expon_dis(mean, size=1, genomeLen=np.inf):
    samples = st.expon.rvs(213.98910256668592,
                           6972.5319847131141, size=size)
    samples = samples * mean / 7106.0
    samples = samples.astype(int)
    samples = np.clip(samples, 1, genomeLen)
    return samples


def draw_mix_gamma_dis(mean, size=1, genomeLen=np.inf):
    """
    lambda
    the actual length should be multiplied  by 1000
    two mixture gamma distribution with parameters
    first gamma: alpha: 6.3693711, rate: 0.53834893
    second gamma: alpha: 1.67638771, rate: 0.22871401
    """
    half = int(size/2.0)
    sample_1 = st.gamma.rvs(6.3693711, 0.53834893, size=half)
    sample = st.gamma.rvs(1.67638771, 0.22871401, size=(size-half))
    sample = np.concatenate((sample, sample_1))
    np.random.shuffle(sample)
    sample = sample*mean/4.39
    sample = sample.astype(int)
    sample = np.clip(sample, 1, genomeLen)
    return sample

#------------------------------------------------------------------------------#


base = ["A", "T", "C", "G"]


def createRandomGenome(fileName="test.genome", length=5000000):
    """
    Creates a random genome of size length and writes it in fileName.genome
    The first line is an integer representing the length of the genome and the
    second line is the entire genome
    """
    with open(fileName, "w") as f:
        f.write("%d\n" % length)
        f.write(''.join(random.choices(base, k=length)))
        f.write('\n')


def reverseComplement(read):
    """
    Returns the reverse complemnt of read
    :param read: a string consisting of only 'A', 'T', 'C', 'G's
    """
    complementDict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complementDict[c] for c in read[::-1]])


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


def generateReads(genomeFileName,
                  readFileName,
                  coverage=25,
                  readLen=10000,
                  pIn=0.00,
                  pDel=0.00,
                  pS=0.05,
                  includeReverse=False,
                  lengthDis=draw_const_dis,
                  ):
    """
    Generate random reads of length readLen with given coverage from the genome
    in the file fileName.genome and stores the reads in fileName.reads. Each
    read is stored in two lines:
    [cn]:readPos:edit string
    readData
    c means reverse [c]omplement; n means [n]ormal
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
    :param includeReverse: whether to include reverse complements
    """
    genome = ''
    genomeLength = 0
    with open(genomeFileName, "r") as f:
        genomeLength = int(f.readline())
        genome = f.readline()[:-1]
    with open(readFileName, "w") as f:
        numTargetBases = genomeLength * coverage
        numCurrentBases = 0
        while numCurrentBases < numTargetBases:
            length = lengthDis(readLen)[0]
            numCurrentBases += length
            (editString, read) = generateRead(genome, length, pIn, pDel, pS)
            reverse = random.choice([True, False]) and includeReverse
            reverseSymbol = 'c' if reverse else 'n'
            f.write(reverseSymbol + ':')
            f.write(editString)
            f.write('\n')
            f.write(reverseComplement(read) if reverse else read)
            f.write('\n')


def generateGenomeFromConfig(section):
    """
    Generates a genome of a specific length from the config section
    """
    length = section.getint("length")
    assert(length is not None)
    genomeFile = section.get("genomeFile")
    assert(genomeFile is not None)
    createRandomGenome(fileName=genomeFile, length=length)


def generateReadsFromConfig(section):
    """
    Generate reads from the config section
    """
    coverage = section.getint("coverage")
    assert(coverage is not None)
    pIn = section.getfloat("pIn")
    assert(pIn is not None)
    pDel = section.getfloat("pDel")
    assert(pDel is not None)
    pS = section.getfloat("pS")
    assert(pS is not None)
    reverseComplement = section.getboolean("reverseComplement")
    assert(reverseComplement is not None)
    readLen = section.getint("readLen")
    assert(readLen is not None)
    lengthDis = section.get("lengthDis")
    assert(lengthDis is not None)
    genomeFile = section.get("genomeFile")
    assert(genomeFile is not None)
    readsFile = section.get("readsFile")
    assert(readsFile is not None)
    lengthDistributions = {
        "const": draw_const_dis,
        "exponential": draw_expon_dis,
        "beta": draw_beta_dis,
        "mixed-gamma": draw_mix_gamma_dis,
    }
    lengthDis = lengthDistributions[lengthDis]
    generateReads(genomeFileName=genomeFile,
                  coverage=coverage,
                  readLen=readLen,
                  pIn=pIn,
                  pDel=pDel,
                  pS=pS,
                  includeReverse=reverseComplement,
                  readFileName=readsFile,
                  lengthDis=lengthDis,
                  )


def generateDataFromConfig(configFileName):
    """
    Generates a random genome and/or random reads based on the configurations 
    in configFileName
    The config file should be in the windows ini format
    An example config file:
    ----------------------------------------------------------------------------
    [genome1]
    genomeLen = 1000000
    genomeFile = test1
    [reads1]
    coverage = 25
    pIn = 0.03
    pDel = 0.03
    pS = 0.04
    reverseComplement = true
    readLen = 10000
    lengthDis = const
    genomeFile = test1
    readsFile = test1.reads
    ----------------------------------------------------------------------------
    Sections should either begin with "genome" or "reads". A section beginning 
    with "genome" will be used to generated a random genome of a certain length 
    and a section beginning with "reads" will be used to generate random reads 
    from a given genome file
    lengthDis specifies the distribution of the read lengths, and should be one 
    of "const", "exponential", "beta", and "mixed-gamma"
    All of the options are REQUIRED.
    """
    config = configparser.ConfigParser()
    config.read(configFileName)
    for section in config.sections():
        if section.startswith("genome"):
            generateGenomeFromConfig(config[section])
        elif section.startswith("reads"):
            generateReadsFromConfig(config[section])


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print("Usage: createData.py initFile")
    else:
        generateDataFromConfig(sys.argv[1])

#------------------------------------------------------------------------------#
# Tests


class TestReverseComplement(unittest.TestCase):

    def test_single(self):
        self.assertEqual(reverseComplement('A'), 'T')
        self.assertEqual(reverseComplement('T'), 'A')
        self.assertEqual(reverseComplement('C'), 'G')
        self.assertEqual(reverseComplement('G'), 'C')

    def test_multiple(self):
        self.assertEqual(reverseComplement('ATTTCAGG'), 'CCTGAAAT')

    def test_long(self):
        length = 1000000
        original = ''.join(random.choices(base, k=length))
        transformed = reverseComplement(original)
        self.assertEqual(reverseComplement(transformed), original)
