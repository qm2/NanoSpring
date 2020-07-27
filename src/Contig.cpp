#include "../include/Contig.h"
#include <omp.h>
#include <algorithm>
#include <stdexcept>
#include <cmath>

ContigGenerator::ContigGenerator(ReadAligner *rA, NanoporeReads &nR, ReadFilter *rF) : rA(rA), nR(nR), rF(rF) {

}

void ContigGenerator::generateContigs() {

    initialize();

    // We start with contig 0, lengthening it at both ends until it is impossible to lengthen it further
    // Then we create a new contig and adds a random read to it
    // If two contigs have overlaps, we merge them

    while (!reads2Contig.empty()) {
        if (!hasActiveContig) {
            // Create a new contig
            read_t r2Add = *reads2Contig.begin();
            reads2Contig.erase(reads2Contig.begin());

            activeContig = new Contig();
            contigs.insert(activeContig);
            activeContig->reads.insert(std::make_pair(0, r2Add));

            readsInContig[r2Add] = std::make_pair(activeContig, 0);
//            std::cout << "Created contig " << contigs.size() - 1 << std::endl;
        }

        hasActiveContig = false;

        {
            size_t l = reads2Contig.size();
            if (l % 20 == 0)
                std::cout << l << " reads left to contig " << contigs.size() << std::endl;
        }

        while (addRelatedReads(*activeContig->reads.begin())) {
            size_t l = reads2Contig.size();
            if (l % 20 == 0)
                std::cout << l << " reads left to contig " << contigs.size() << std::endl;
        }

        while (true) {
            auto lastRead = activeContig->reads.end();
            lastRead--;
            if (!addRelatedReads(*lastRead))
                break;
            size_t l = reads2Contig.size();
            if (l % 20 == 0)
                std::cout << l << " reads left to contig " << contigs.size() << std::endl;
        }
    }
}

void ContigGenerator::initialize() {
    hasActiveContig = false;
    readsInContig.clear();
    reads2Contig.clear();
    for (auto c: contigs)
        delete c;
    contigs.clear();
    std::vector<read_t> temp;
    temp.reserve(nR.numReads);
    for (size_t i = 0; i < nR.numReads; ++i)
        temp.push_back(i);
    std::random_shuffle(temp.begin(), temp.end());
    reads2Contig.insert(temp.begin(), temp.end());
}

bool ContigGenerator::addRelatedReads(const std::pair<long, read_t> r) {
//            std::cout << "Adding related reads" << std::endl;
    std::vector<size_t> results;
    rF->getFilteredReads(r.second, results);
    bool addedRead = false;
    size_t const resultLen = results.size();
    bool merged = false;
    long relPosInMerge;
    Contig *contig2MergeWith;

    omp_lock_t lock;
    omp_init_lock(&lock);
#pragma omp parallel for shared(lock, results, addedRead, merged, contig2MergeWith, relPosInMerge)
    for (size_t i = 0; i < resultLen; i++) {
        size_t it = results[i];
        long relPos;
        omp_set_lock(&lock);
        auto findRead = readsInContig.find(it);
        if (findRead != readsInContig.end()) {
            std::pair<Contig *, long> matchedRead = findRead->second;
            if (!merged && matchedRead.first != activeContig &&
                rA->align(*nR.readData[r.second], *nR.readData[it], relPos)) {
                merged = true;
                contig2MergeWith = matchedRead.first;
                relPosInMerge = matchedRead.second - relPos - r.first;
            }
            omp_unset_lock(&lock);
            continue;
        } else {
            omp_unset_lock(&lock);
        }

        if (rA->align(*nR.readData[r.second], *nR.readData[it], relPos)) {
            long pos2Add = r.first + relPos;
            {
                omp_set_lock(&lock);
                activeContig->reads.insert(std::make_pair(pos2Add, it));
                reads2Contig.erase(it);
                readsInContig[it] = std::make_pair(activeContig, pos2Add);
                addedRead = true;
                omp_unset_lock(&lock);
            }
        }
    }
    omp_destroy_lock(&lock);

    if (merged) {
        // Merges the second contig into the first one
        if (contig2MergeWith->reads.size() < activeContig->reads.size())
            mergeContigs(activeContig, contig2MergeWith, -relPosInMerge);
        else
            mergeContigs(contig2MergeWith, activeContig, relPosInMerge);
    }
    return addedRead;
};

void ContigGenerator::mergeContigs(Contig *c1, Contig *c2, long pos) {
//                    std::cout << "Merging" << std::endl;
    long smallestPos = c1->reads.begin()->first;
    auto it = c1->reads.end();
    it--;
    long biggestPos = it->first;
    for (const std::pair<long, size_t> &r : c2->reads) {
        long pos2Add = r.first + pos;
        if (pos2Add > biggestPos)
            hasActiveContig = true;
        else if (pos2Add < smallestPos)
            hasActiveContig = true;
        c1->reads.insert(std::make_pair(pos2Add, r.second));
        // Update the map from reads to contigs
        readsInContig[r.second] = std::make_pair(c1, pos2Add);
//                        std::cout << r.second << " ";
    }
//                    std::cout << std::endl;
    activeContig = c1;
    hasActiveContig = true;
    contigs.erase(c2);
    delete c2;
}

std::ostream &operator<<(std::ostream &out, const ContigGenerator &o) {
    size_t l = o.contigs.size();
    out << l << " contigs are generated" << std::endl;
    size_t i = 0;
    for (auto contig : o.contigs) {
        size_t contigSize = contig->reads.size();
        out << contigSize << " reads in contig " << i++ << std::endl;
        double averageShift = 0;
        size_t j = 0;

        const long genomeLen = 100000;
//        const long genomeLen = 5000000;
        for (const auto r : contig->reads) {
            long dif = r.first - (long) o.nR.readPos[r.second];
            dif += 2 * genomeLen;
            dif %= genomeLen;
            if (dif > genomeLen / 2)
                dif -= genomeLen;
            averageShift += dif;
//            out << r.first << " " << (long) o.nR.readPos[r.second] << std::endl;
        }

        averageShift /= contigSize;

        long difErrs[contigSize];
        for (const auto r : contig->reads) {
            long dif = r.first - (long) o.nR.readPos[r.second];
            dif += 2 * genomeLen;
            dif %= genomeLen;
            if (dif > genomeLen / 2)
                dif -= genomeLen;
            difErrs[j++] = dif - averageShift;
//            out << r.first << " " << (long) o.nR.readPos[r.second]
//                << " " << dif - averageShift << " " << r.second << std::endl;
        }

        std::sort(difErrs, difErrs + contigSize);

        out << "Average shift is " << averageShift << std::endl;
        out << "Shift Error median is " << difErrs[contigSize / 2] << std::endl;
        out << "Shift Error 25% is " << difErrs[contigSize / 4] << std::endl;
        out << "Shift Error 75% is " << difErrs[(contigSize * 3) / 4] << std::endl;
    }
    return out;
}

ContigGenerator::~ContigGenerator() {
    for (auto c: contigs)
        delete c;
}