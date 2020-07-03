import unittest, numpy as np

class NanoRead:
    """
    An instance of class NanoRead is one read sequence
    """

    def __init__(self, line1, line2):
        """
        Initializes a read from two lines in the .reads file
        """
        line1Split = line1.split(':')
        self.pos = int(line1Split[0])
        self.editString = line1Split[1]
        self.data = line2

    def getMinHash(self, k, n):
        """
        Returns a numpy array of size n which is a MinHash sketch of the read
        """
        try:
            if self.k == k:
                if self.n == n:
                    return self.minHash
            self.k = k
            self.n = n
            self.minHash = NanoRead.__calcMinHash(self.k, self.n)
            return self.minHash
        except AttributeError:
            self.k = k
            self.n = n
            self.minHash = NanoRead.__calcMinHash(self.k, self.n)
            return self.minHash

    def __calcMinHash(self, k, n):
        """
        Returns a numpy array of size n which is a MinHash sketch of the read
        """
        l = len(self.data)
        numKMers = l - k + 1
        hashes = np.zeros((numKMers, n), dtype=(np.uint32))
        for i in range(numKMers):
            hashes[i, :] = NanoRead.calcHash(self.data[i:i + k], n)

        return np.amin(hashes, axis=0)

    @staticmethod
    def calcHash(kMer, n):
        """
        Returns a numpy array of n different hashes of kMer. After the first
        hash is calculated, the other hashes are just bit rotations + xor
        """
        hashes = np.zeros((n,), dtype=(np.uint32))
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


class TestNanoRead(unittest.TestCase):

    def testInit(self):
        pos = 339
        editString = 'u4dAsTG'
        data = 'ATTTGC'
        n = NanoRead(str(pos) + ':' + editString, data)
        self.assertEqual(n.pos, pos)
        self.assertEqual(n.editString, editString)
        self.assertEqual(n.data, data)

    def testMinHash(self):
        pos = 339
        editString = 'u4dAsTG'
        data1 = 'ATTTGCGCCCTACCGT'
        data2 = 'ATCTGCACCCTACCCC'
        nR1 = NanoRead(str(pos) + ':' + editString, data1)
        nR2 = NanoRead(str(pos) + ':' + editString, data2)
        sketch1 = nR1.getMinHash(3, 4)
        sketch2 = nR2.getMinHash(3, 4)
        print(sketch1 == sketch2)

    def testCalcHash(self):
        kMer = 'ACTTGGG'
        print(NanoRead.calcHash(kMer, 8))
        print(NanoRead.calcHash(kMer, 8))


if __name__ == '__main__':
    unittest.main()
