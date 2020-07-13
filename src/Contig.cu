//
// Created by yifan on 7/10/20.
//

#include "../include/Contig.cuh"

ContigGenerator::ContigGenerator(ReadAligner *rA, NanoporeReads &nR, ReadFilter *rF) : rA(rA), nR(nR), rF(rF) {

}

void ContigGenerator::generateContigs() {
    // We maintain a list of all reads that have not been treated and a map from reads to contigs
    std::unordered_set<read_t> reads2Contig;
    for (size_t i = 0; i < nR.numReads; ++i)
        reads2Contig.insert(i);
    std::map<read_t, size_t> readsInContig;



    // We start with contig 0, lengthening it at both ends until it is impossible to lengthen it further
    // Then we create a new contig and adds a random read to it
    size_t currentContigIndex = (size_t) -1;

    // We create a contig and add the first read to it.
    auto initialize = [&]() {
        read_t r = *reads2Contig.begin();
        reads2Contig.erase(reads2Contig.begin());

        contigs.push_back(Contig());
        contigs[++currentContigIndex].reads.insert(std::make_pair(0, r));

        readsInContig.insert(std::make_pair(r, currentContigIndex));
        std::cout << "Created contig " << currentContigIndex << std::endl;
    };

    while (!reads2Contig.empty()) {
        initialize();
        Contig &currentContig = contigs[currentContigIndex];

        auto addRelatedReads = [&](std::pair<long, read_t> r) {
            std::vector<size_t> results;
            rF->getFilteredReads(r.second, results);
            bool addedRead = false;
            for (const auto &it : results) {
                long relPos;
                if (reads2Contig.find(it) == reads2Contig.end()) {
                    continue;
                }
                if (rA->align(*nR.readData[r.second], *nR.readData[it], relPos)) {
                    currentContig.reads.insert(std::make_pair(r.first + relPos, it));
                    reads2Contig.erase(it);
                    addedRead = true;
                }
            }
            return addedRead;
        };


        while (addRelatedReads(*currentContig.reads.begin())) {
            size_t l = reads2Contig.size();
            std::cout << l << " reads left to contig" << std::endl;
        }

        while (true) {
            auto lastRead = currentContig.reads.end();
            lastRead--;
            if (!addRelatedReads(*lastRead))
                break;
            size_t l = reads2Contig.size();
            std::cout << l << " reads left to contig" << std::endl;
        }
    }
//    while (!reads2Contig.empty()) {
//
//    }
}

std::ostream &operator<<(std::ostream &out, const ContigGenerator &o) {
    size_t l = o.contigs.size();
    out << l << " contigs are generated" << std::endl;
    for (size_t i = 0; i < l; ++i) {
        out << o.contigs[i].reads.size() << " reads in contig " << i << std::endl;
    }
    return out;
}
