//
// Created by yifan on 7/10/20.
//

#include "../include/Contig.cuh"

ContigGenerator::ContigGenerator(ReadAligner *rA, NanoporeReads &nR, ReadFilter *rF) : rA(rA), nR(nR), rF(rF) {

}

void ContigGenerator::generateContigs() {
    // We maintain a list of all reads that have not been treated and a map from reads to contigs
    std::unordered_set<read_t> reads2Contig;
    {
        std::vector<read_t> temp;
        temp.reserve(nR.numReads);
        for (size_t i = 0; i < nR.numReads; ++i)
            temp.push_back(i);
        std::random_shuffle(temp.begin(), temp.end());
        reads2Contig.insert(temp.begin(), temp.end());
    }
    std::map<read_t, std::pair<size_t, long>> readsInContig;



    // We start with contig 0, lengthening it at both ends until it is impossible to lengthen it further
    // Then we create a new contig and adds a random read to it
    size_t currentContigIndex = (size_t) -1;

    while (!reads2Contig.empty()) {
        {
            read_t r2Add = *reads2Contig.begin();
            std::vector<size_t> results;
            rF->getFilteredReads(r2Add, results);
            bool successfullyAdded = false;
            for (const auto &it : results) {
                long relPos;
                // If "it" has already been added to a contig
                if (reads2Contig.find(it) == reads2Contig.end()) {
                    // If there is alignment
                    if (rA->align(*nR.readData[r2Add], *nR.readData[it], relPos)) {
                        size_t contigId = readsInContig[it].first;
                        Contig &contig2Add = contigs[contigId];
                        long pos2Add = readsInContig[it].second + relPos;
                        contig2Add.reads.insert(std::make_pair(pos2Add, r2Add));
                        reads2Contig.erase(r2Add);
                        readsInContig[r2Add] = std::make_pair(contigId, pos2Add);
                        successfullyAdded = true;
                        break;
                    }
                }
            }
            {
                size_t l = reads2Contig.size();
                if (l % 20 == 0)
                    std::cout << l << " reads left to contig[Adding to existing contigs]" << std::endl;
            }
            if (successfullyAdded)
                continue;
            // Otherwise create a new contig
            reads2Contig.erase(reads2Contig.begin());

            contigs.push_back(Contig());
            contigs[++currentContigIndex].reads.insert(std::make_pair(0, r2Add));

            readsInContig[r2Add] = std::make_pair(currentContigIndex, 0);
            std::cout << "Created contig " << currentContigIndex << std::endl;
        }

        Contig &currentContig = contigs[currentContigIndex];

        auto addRelatedReads = [&](std::pair<long, read_t> r) {
            std::vector<size_t> results;
            rF->getFilteredReads(r.second, results);
            bool addedRead = false;
            size_t resultLen = results.size();
            omp_lock_t lock;
            omp_init_lock(&lock);
#pragma omp parallel for
            for (size_t i = 0; i < resultLen; i++) {
                size_t it = results[i];
                long relPos;
                omp_set_lock(&lock);
                if (reads2Contig.find(it) == reads2Contig.end()) {
                    omp_unset_lock(&lock);
                    continue;
                }
                omp_unset_lock(&lock);
                if (rA->align(*nR.readData[r.second], *nR.readData[it], relPos)) {
                    long pos2Add = r.first + relPos;
                    {
                        omp_set_lock(&lock);
                        currentContig.reads.insert(std::make_pair(pos2Add, it));
                        reads2Contig.erase(it);
                        readsInContig[it] = std::make_pair(currentContigIndex, pos2Add);
                        addedRead = true;
                        omp_unset_lock(&lock);
                    }
                }
            }
            omp_destroy_lock(&lock);
            return addedRead;
        };


        while (addRelatedReads(*currentContig.reads.begin())) {
            size_t l = reads2Contig.size();
            if (l % 10 == 0)
                std::cout << l << " reads left to contig" << std::endl;
        }

        while (true) {
            auto lastRead = currentContig.reads.end();
            lastRead--;
            if (!addRelatedReads(*lastRead))
                break;
            size_t l = reads2Contig.size();
            if (l % 10 == 0)
                std::cout << l << " reads left to contig" << std::endl;
        }
    }
}

std::ostream &operator<<(std::ostream &out, const ContigGenerator &o) {
    size_t l = o.contigs.size();
    out << l << " contigs are generated" << std::endl;
    for (size_t i = 0; i < l; ++i) {
        size_t contigSize = o.contigs[i].reads.size();
        out << contigSize << " reads in contig " << i << std::endl;
        double averageShift = 0;
        size_t j = 0;

        for (const auto r : o.contigs[i].reads) {
            long dif = r.first - (long) o.nR.readPos[r.second];
            const long genomeLen = 5000000;
            dif += 2 * genomeLen;
            dif %= genomeLen;
            if (dif > genomeLen / 2)
                dif -= genomeLen;
            averageShift += dif;
//            out << r.first << " " << (long) o.nR.readPos[r.second] << std::endl;
        }

        averageShift /= contigSize;

        long difErrs[contigSize];
        for (const auto r : o.contigs[i].reads) {
            long dif = r.first - (long) o.nR.readPos[r.second];
            const long genomeLen = 5000000;
            dif += 2 * genomeLen;
            dif %= genomeLen;
            if (dif > genomeLen / 2)
                dif -= genomeLen;
            difErrs[j++] = dif - averageShift;
//            out << r.first << " " << (long) o.nR.readPos[r.second] << std::endl;
        }

        std::sort(difErrs, difErrs + contigSize);

        out << "Average shift is " << averageShift << std::endl;
        out << "Shift Error median is " << difErrs[contigSize / 2] << std::endl;
        out << "Shift Error 25% is " << difErrs[contigSize / 4] << std::endl;
        out << "Shift Error 75% is " << difErrs[(contigSize * 3) / 4] << std::endl;
    }
    return out;
}
