import unittest
import numpy as cp

class NanoRead:
    """
    An instance of class NanoRead is one read sequence
    """

    def __init__(self, index, line1, line2):
        """
        Initializes a read from two lines in the .reads file
        """
        line1Split = line1.split(':')
        self.pos = int(line1Split[0])
        self.editString = line1Split[1]
        self.index = index
        self.data = line2

    def __str__(self):
        return "%d %d:%.10s %.10s" % (self.index, self.pos, self.editString, self.data)

    def getMinHash(self, k, n):
        """
        Returns a numpy array of size n which is a MinHash sketch of the read.
        For speed reasons don't use this!
        """
        if (self.index % 100 == 0):
            print("Building minHash of read " + str(self))
        try:
            if self.k == k:
                if self.n == n:
                    return self.minHash
            self.k = k
            self.n = n
            self.minHash = self.__calcMinHash(self.k, self.n)
            return self.minHash
        except AttributeError:
            self.k = k
            self.n = n
            self.minHash = self.__calcMinHash(self.k, self.n)
            return self.minHash

    def __calcMinHash(self, k, n):
        """
        Returns a numpy array of size n which is a MinHash sketch of the read
        """
        l = len(self.data)
        numKMers = l - k + 1
        hashes = cp.zeros((numKMers, n), dtype=(cp.uint32))
        for i in range(numKMers):
            hashes[i, :] = NanoRead.calcHash(self.data[i:i + k], n)

        return cp.amin(hashes, axis=0)

    @staticmethod
    def calcHash(kMer, n):
        """
        Returns a numpy array of n different hashes of kMer. After the first
        hash is calculated, the other hashes are just bit rotations + xor
        """
        hashes = cp.zeros((n,), dtype=(cp.uint32))
        hashes[0] = hash(kMer)
        INT_BITS = 32
        # Number of bits to rotate left
        d = 5
        for i in range(1, n):
            hashes[i] = hashes[(i - 1)] << d | hashes[(i - 1)] >> INT_BITS - d
            hashes[i] ^= 0xABCD0123

        return hashes

class NanoReads:
    """
    Stores all reads from a genome
    """

    def __init__(self, fileName):
        """
        Initializes the data from fileName.reads
        """
        self.reads = []
        i = 0
        with open(fileName + ".reads") as f:
            while (True):
                l = f.readline()
                if (l == ''):
                    break
                self.reads.append(NanoRead(i, l[:-1], f.readline()[:-1]))
                i += 1
        self.sortedReads = self.reads[:]
        self.sortedReads.sort(key = lambda nR : nR.pos)
        self.readLen = len(self.reads[0].data)
        self.numReads = len(self.reads)

    def buildHashTables(self, k, n):
        """
        Builds n Hash Tables with the sketches obtained by MinHash on k-mers
        """
        self.sketches = cp.array([nR.getMinHash(k, n) for nR in self.reads])
        self.hashTables = []
        for i in range(n):
            self.hashTables.append({})
            for j in range(self.numReads):
                self.hashTables[i][self.sketches[j, i]] = self.reads[j]
            


class TestNanoRead(unittest.TestCase):

    def testInit(self):
        pos = 339
        editString = 'u4dAsTG'
        data = 'ATTTGC'
        n = NanoRead(0, str(pos) + ':' + editString, data)
        self.assertEqual(n.pos, pos)
        self.assertEqual(n.editString, editString)
        self.assertEqual(n.data, data)

    def testMinHash(self):
        pos = 339
        editString = 'u4dAsTG'
        data1 = 'ATTTGCGCCCTACCGT'
        data2 = 'ATCTGCACCCTACCCC'
        nR1 = NanoRead(0, str(pos) + ':' + editString, data1)
        nR2 = NanoRead(1, str(pos) + ':' + editString, data2)
        sketch1 = nR1.getMinHash(3, 4)
        sketch2 = nR2.getMinHash(3, 4)
        print(sketch1 == sketch2)

    def testCalcHash(self):
        kMer = 'ACTTGGG'
        print(NanoRead.calcHash(kMer, 8))
        print(NanoRead.calcHash(kMer, 8))

    def testNanoReads(self):
        fileName = "test1"
        nanoReads = NanoReads(fileName)
        #print('\n'.join(map(str, nanoReads.reads)))
        #print('\n'.join(map(str, nanoReads.sortedReads)))
        #print(nanoReads.readLen)
        print("Building Hash Tables");
        nanoReads.buildHashTables(9, 8)

if __name__ == '__main__':
    unittest.main()
