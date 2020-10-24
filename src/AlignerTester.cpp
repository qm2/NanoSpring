#include "AlignerTester.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <sstream>

void AlignerTester::generateData(size_t readLen, ssize_t offset, size_t num,
                                 double pIns, double pDel, double pS) {
    readsA.clear();
    readsB.clear();

    const std::string BASES = "ATCG";

    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_int_distribution<> baseDis(0, BASES.length() - 1);

    for (size_t i = 0; i < num; ++i) {
        // First we generate the original sequence. THe original sequence does
        // not end in \0!
        size_t origSeqLen = readLen * 2;
        char origSeq[origSeqLen];
        for (size_t j = 0; j < origSeqLen; ++j)
            origSeq[j] = BASES[baseDis(generator)];

        // The function for generating the read and adding it to reads
        auto generateRead = [&origSeq, origSeqLen, readLen, pIns, pDel, pS,
                             &generator, &baseDis,
                             &BASES](ssize_t offset,
                                     std::vector<std::string> &reads) {
            char *origSeqEnd = origSeq + origSeqLen;
            char *currPos =
                offset >= 0 ? origSeq + offset : origSeqEnd + offset;
            std::string read;
            size_t k = 0;
            std::uniform_real_distribution<> realDis(0, 1);
            while (k < readLen) {
                double r = realDis(generator);
                if (r < pIns) {
                    // Insert a random base
                    read.push_back(BASES[baseDis(generator)]);
                    k++;
                } else if (r >= pIns && r < pIns + pDel) {
                    // Delete a base
                    currPos++;
                    if (currPos == origSeqEnd)
                        currPos = origSeq;
                } else if (r >= pIns + pDel && r < pIns + pDel + pS) {
                    // Substitution error
                    char oldBase = *currPos;
                    std::string basesToChoose = BASES;
                    basesToChoose.erase(std::remove(basesToChoose.begin(),
                                                    basesToChoose.end(),
                                                    oldBase),
                                        basesToChoose.end());
                    std::uniform_int_distribution<> subBaseDis(
                        0, basesToChoose.length() - 1);
                    read.push_back(basesToChoose[subBaseDis(generator)]);
                    k++;
                    currPos++;
                    if (currPos == origSeqEnd)
                        currPos = origSeq;

                } else {
                    // No error
                    read.push_back(*currPos);
                    k++;
                    currPos++;
                    if (currPos == origSeqEnd)
                        currPos = origSeq;
                }
            }
            reads.push_back(read);
            //            std::cout << read << std::endl;
        };
        generateRead(0, readsA);
        generateRead(offset, readsB);
    }
}

void AlignerTester::profile(StringAligner<const char *> *aligner,
                            double &duration, double &successRate,
                            double &avgBeginOffset, double &avgEndOffset,
                            double &avgEditDis) {
    size_t numSuccess = 0;
    size_t totalEditDis = 0;
    ssize_t totalBeginOffset = 0;
    ssize_t totalEndOffset = 0;
    auto start = std::chrono::high_resolution_clock::now();

    size_t numReads = readsA.size();
#pragma omp parallel for reduction(+:numSuccess, totalEditDis, totalBeginOffset, totalEndOffset)
    for (size_t i = 0; i < numReads; ++i) {
        std::vector<Edit> editScript;
        size_t editDis;
        ssize_t beginOffset, endOffset;
        const char *Abegin = readsA[i].c_str();
        const char *Aend = Abegin + readsA[i].length();
        const char *Bbegin = readsB[i].c_str();
        const char *Bend = Bbegin + readsB[i].length();
        bool success =
            aligner->align(Abegin, Aend, Bbegin, Bend, 0, beginOffset,
                           endOffset, editScript, editDis);
        if (success) {
            numSuccess++;
            totalEditDis += editDis;
            totalBeginOffset += beginOffset;
            totalEndOffset += endOffset;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000000000.0 / (double)numReads;
    successRate = numSuccess / (double)numReads;
    avgEditDis = totalEditDis / (double)numSuccess;
    avgBeginOffset = totalBeginOffset / (double)numSuccess;
    avgEndOffset = totalEndOffset / (double)numSuccess;
}

void AlignerTester::profile(StringAligner<const char *> *aligner,
                            double &duration, double &successRate,
                            double &avgEditDis) {
    size_t numSuccess = 0;
    size_t totalEditDis = 0;
    auto start = std::chrono::high_resolution_clock::now();

    size_t numReads = readsA.size();
    for (size_t i = 0; i < numReads; ++i) {
        std::vector<Edit> editScript;
        size_t editDis;
        const char *Abegin = readsA[i].c_str();
        const char *Aend = Abegin + readsA[i].length();
        const char *Bbegin = readsB[i].c_str();
        const char *Bend = Bbegin + readsB[i].length();
        ssize_t beginOffset, endOffset;
        bool success =
            aligner->align(Abegin, Aend, Bbegin, Bend, 0, beginOffset,
                           endOffset, editScript, editDis);
        if (success) {
            numSuccess++;
            totalEditDis += editDis;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000000000.0 / (double)numReads;
    successRate = numSuccess / (double)numReads;
    avgEditDis = totalEditDis / (double)numSuccess;
}

bool AlignerTester::validate(StringAligner<const char *> *aligner) {
    size_t numReads = readsA.size();
    bool result = true;
#pragma omp parallel
    {
#pragma omp for
        for (size_t i = 0; i < numReads; ++i) {
            if (!result)
                continue;
            std::vector<Edit> editScript;
            size_t editDis;
            ssize_t beginOffset, endOffset;
            // std::cout << readsA[i] << '\n';
            // std::cout << readsB[i] << '\n';

            const char *Abegin = readsA[i].c_str();
            const char *Aend = Abegin + readsA[i].length();
            const char *Bbegin = readsB[i].c_str();
            const char *Bend = Bbegin + readsB[i].length();
            bool success =
                aligner->align(Abegin, Aend, Bbegin, Bend, 0, beginOffset,
                               endOffset, editScript, editDis);
            if (!success) {
                result = false;
#pragma omp cancel for
                continue;
            }
            std::string resultAfterEdit;
            // std::cout << "beginOffset" << beginOffset << std::endl;
            // std::cout << "endOffset" << endOffset << std::endl;
            // std::cout << aligner->name << std::endl;
            std::string origString = readsA[i].substr(
                beginOffset > 0 ? beginOffset : 0,
                readsA[i].length() - (beginOffset > 0 ? beginOffset : 0) +
                    (endOffset > 0 ? 0 : endOffset));
            std::string targetString = readsB[i].substr(
                beginOffset > 0 ? 0 : -beginOffset,
                readsB[i].length() - (beginOffset > 0 ? 0 : -beginOffset) -
                    (endOffset > 0 ? endOffset : 0));

            Edits::applyEdits(
                origString.begin(), editScript,
                std::inserter(resultAfterEdit, resultAfterEdit.end()));
            if (resultAfterEdit.compare(targetString)) {
                std::cout << "Begin Offset " << beginOffset << " EndOffset "
                          << endOffset << '\n';
                std::cout << readsA[i] << std::endl;
                std::cout << readsB[i].length() << readsB[i] << std::endl;
                for (auto e : editScript)
                    std::cout << e;
                std::cout << std::endl;
                std::cout << "result is " << result << std::endl;
                std::cout << resultAfterEdit.length() << resultAfterEdit
                          << std::endl;
                result = false;
#pragma omp cancel for
            }
        }
    }
    return result;
}