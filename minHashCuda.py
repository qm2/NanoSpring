import unittest
import numpy as cp

class NanoReads:
    """
    Stores all reads from a genome
    """

    def __init__(self, fileName):
        """
        Initializes the data from fileName.reads
        """
        self.data = []
        self.pos = []
        self.editStrings = []
        i = 0
        with open(fileName + ".reads") as f:
            while (True):
                l = f.readline()
                if (l == ''):
                    break
                lSplit = l.split(':')
                self.pos.append(int(lSplit[0]))
                self.editStrings.append(lSplit[1])
                i += 1
        self.sortedReads = self.reads[:]
        self.sortedReads.sort(key = lambda nR : nR.pos)
        self.readLen = len(self.reads[0].data)
        self.numReads = len(self.reads)


class TestNanoRead(unittest.TestCase):

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
