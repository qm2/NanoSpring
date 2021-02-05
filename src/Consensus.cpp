#include "Consensus.h"
#include "DirectoryUtils.h"
#include "bsc_helper.h"
#include <algorithm>
#include <arpa/inet.h>
#include <boost/filesystem.hpp>
#include <csignal>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <set>
#include "minimap.h"

void Consensus::generateAndWriteConsensus() {
    initialize();

    std::vector<std::vector<read_t>> numReadsInContig(numThr);
    std::vector<std::vector<read_t>> loneReads(numThr);
    ConsensusGraph *cG = nullptr;
    std::vector<CountStats> count_stats(numThr);
#pragma omp parallel private(cG)
    {
        auto tid = omp_get_thread_num();
        std::ofstream logfile("logfile"+std::to_string(tid), std::ofstream::out);
        std::string filePrefix = tempDir + tempFileName + std::to_string(tid);
        ConsensusGraphWriter cgw(filePrefix);
        read_t firstUnaddedRead = 0; 
        int contigId = 0;
        // guarantee that all reads < firstUnaddedRead have been picked
        while ((cG = createGraph(firstUnaddedRead))) {
            logfile<<"Thread: " << omp_get_thread_num() << ", Contig: " << contigId << ", First read number "<<cG->readsInGraph.begin()->first<<"\n";
            ssize_t initialStartPos = cG->startPos; // simply 0
            ssize_t initialEndPos = cG->endPos; // 
            const ssize_t len = initialEndPos - initialStartPos;
            // TODO: Experiment with this
            size_t offset = rD->avgReadLen / 4;

            ssize_t curPos = cG->startPos;

            while (true) {
                std::cout << "right\n";
                addRelatedReads(cG, curPos, len, count_stats[tid], logfile, contigId);
                cG->printStatus();
                curPos += offset;
                // std::cout << "curPos " << curPos << " len " << len << " endPos "
                //          << cG->endPos << '\n';
                if (curPos + len > cG->endPos)
                    break;
            }

            curPos = initialStartPos - offset;
            while (true) {
                std::cout << "left\n";
                if (curPos < cG->startPos)
                    break;
                cG->printStatus();
                addRelatedReads(cG, curPos, len, count_stats[tid], logfile, contigId);
                curPos -= offset;
            }
            cG->writeMainPath(cgw);
            cG->writeReads(cgw);
            numReadsInContig[tid].push_back(cG->getNumReads());
            if (cG->getNumReads() == 1)
                loneReads[tid].push_back(cG->readsInGraph.begin()->first);
            delete cG;
            contigId++;
        }
    } // pragma omp parallel

    // now perform last step, combining files from threads and writing metadata
    finishWriteConsensus(numReadsInContig);

    // print stats about lone reads
    read_t totalNumLoneReads = 0;
    std::cout << "LoneReads";
    for (auto &v: loneReads) {
        totalNumLoneReads += v.size();
        for (auto &r: v)
           std::cout << std::dec << ":" << r; 
    }

    // compute total count stats by adding for all threads
    CountStats summary;
    for (auto &c : count_stats)
        summary = summary + c;

    std::cout << "\n";
    std::cout << "#LoneReads = " << totalNumLoneReads << "\n";
    std::cout << "MinHash passed " << summary.countMinHash << std::dec << " reads\n";
    std::cout << "MinHash passed & not already in graph " << summary.countMinHashNotInGraph << std::dec << " reads\n";
    std::cout << "Merge Sort passed " << summary.countMergeSort << std::dec << " reads\n";
    std::cout << "Aligner passed " << summary.countAligner << std::dec << " reads\n";
}

void Consensus::addRelatedReads(ConsensusGraph *cG, ssize_t curPos, int len, CountStats &cs, std::ofstream &logfile, int contigId) {
    // Find reads likely to have overlaps
    ssize_t offsetInMainPath = curPos - cG->startPos;
    if (len == 0 || offsetInMainPath < 0 ||
        offsetInMainPath >= (ssize_t)cG->mainPath.path.size())
        return;
    auto stringBegin = cG->mainPath.path.begin() + offsetInMainPath;
    auto stringEnd =
        (ssize_t)cG->mainPath.path.size() >= offsetInMainPath + (ssize_t)len
            ? stringBegin + len
            : cG->mainPath.path.end();
    // std::cout << curPos << "\n";
    const std::string originalString(stringBegin, stringEnd);
    std::string reverseComplementString;
    ReadData::toReverseComplement(
        stringBegin, stringEnd,
        std::inserter(reverseComplementString, reverseComplementString.end()));
    bool all[] = {false, true};
    for (bool reverseComplement : all) {
        std::vector<read_t> results;
        rF->getFilteredReads(reverseComplement ? reverseComplementString
                                               : originalString,
                             results);
        cs.countMinHash += results.size();
        // std::cout << "Found " << results.size() << " reads\n";
        // for (read_t r : results)
        // std::cout << r << " ";
        // std::cout << '\n';

        // Try to add them one by one
        for (const auto r : results) {
            // std::cout << "Working on read " << r << '\n';
            if (inGraph[r])
                continue;

            cs.countMinHashNotInGraph++;
            logfile<<"Thread: " << omp_get_thread_num() << ", Contig: " << contigId << ", Read passed MinHash "<<r<<"\n";
            ssize_t relPos;
            std::string readStr;
            if (reverseComplement)
                ReadData::toReverseComplement(
                    rD->getRead(r).begin(), rD->getRead(r).end(),
                    std::inserter(readStr, readStr.end()));
            else
                readStr = rD->getRead(r);

            if (!rA->align(originalString, readStr, relPos)) {
                logfile<<"Thread: " << omp_get_thread_num() << ", Contig: " << contigId << ", Read failed Sort-Merge "<<r<<"\n";
                continue;
            }
            logfile<<"Thread: " << omp_get_thread_num() << ", Contig: " << contigId << ", Read passed Sort-Merge "<<r<<"\n";
            cs.countMergeSort++;
            ssize_t pos = curPos + relPos;
            std::vector<Edit> editScript;
            ssize_t beginOffset, endOffset;
            if (!readStatusLock[r%numLocks].try_lock()) {
                // we only try_lock here since missing a read
                // is not a major issue and lock contention should be a rare event anyway.
                // On the other hand, waiting for another thread to finish alignment is not
                // worthwhile.
                continue;
            } else {
                // check again that read is available (variables flushed after lock is set)
                if (inGraph[r]) {
                    // read already taken, continue with next read
                    readStatusLock[r%numLocks].unlock();
                    continue;
                }
                if (!cG->addRead(readStr, pos, editScript, beginOffset,
                                 endOffset)) {
                    // read doesn't align, continue with next read
                    logfile<<"Thread: " << omp_get_thread_num() << ", Contig: " << contigId << ", Read failed aligner "<<r<<"\n";
                    readStatusLock[r%numLocks].unlock();
                    continue;
                } else {
                    // read added to graph
                    inGraph[r] = true;
                    readStatusLock[r%numLocks].unlock();
                    cs.countAligner++;                    
                }
            }
#ifdef CHECKS
            {
                // Check editScript applied to originalString is readStr
                std::string origString = std::string(
                    cG->mainPath.path.begin() +
                        (beginOffset > 0 ? beginOffset : 0),
                    cG->mainPath.path.end() + (endOffset > 0 ? 0 : endOffset));
                std::string targetString = readStr.substr(
                    beginOffset > 0 ? 0 : -beginOffset,
                    readStr.length() - (beginOffset > 0 ? 0 : -beginOffset) -
                        (endOffset > 0 ? endOffset : 0));
                std::string resultAfterEdit;
                Edits::applyEdits(
                    origString.begin(), editScript,
                    std::inserter(resultAfterEdit, resultAfterEdit.end()));
                if (resultAfterEdit.compare(targetString)) {
                    std::cout << "beginOffset " << beginOffset << " endOffset "
                              << endOffset << std::endl;
                    std::cout << "origString\n" << origString << std::endl;
                    std::cout << "targetString\n" << targetString << std::endl;
                    std::cout << "resultAfterEdit\n"
                              << resultAfterEdit << std::endl;
                    std::cout << "editScript\n";
                    for (auto e : editScript)
                        std::cout << e;
                    std::cout << std::endl;
                    std::cout << "mainPath\n"
                              << std::string(cG->mainPath.path.begin(),
                                             cG->mainPath.path.end())
                              << std::endl;
                    std::cout << "pos " << pos << " endPos " << cG->endPos
                              << " offsetGuess "
                              << cG->mainPath.path.size() - cG->endPos + pos
                              << std::endl;
                }
                assert(!resultAfterEdit.compare(targetString));
            }
#endif
            cG->updateGraph(readStr, editScript, beginOffset, endOffset, r, pos,
                            reverseComplement);
#ifdef CHECKS
            assert(checkRead(cG, r));
            assert(cG->checkNoCycle());
#endif
            cG->calculateMainPathGreedy();
#ifdef CHECKS
            std::cout << "Added read " << r << " first unadded read "
                      << firstUnaddedRead << std::endl;
            assert(checkRead(cG, r));
            assert(cG->checkNoCycle());
#endif
        }
    }
}

bool Consensus::checkRead(ConsensusGraph *cG, read_t read) {
    std::string result;
    bool temp = cG->getRead(read, std::inserter(result, result.end()));
    assert(temp);
    if (!temp)
        return false;
    if (!cG->readsInGraph.at(read).reverseComplement) {
        if (result != rD->getRead(read)) {
            std::cout << "readInGraph:\n" << result << "\n";
            std::cout << "actualRead:\n" << rD->getRead(read) << "\n";
        }
        return result == rD->getRead(read);
    }
    std::string reverseComplement;
    ReadData::toReverseComplement(
        result.begin(), result.end(),
        std::inserter(reverseComplement, reverseComplement.end()));
    if (reverseComplement != rD->getRead(read)) {
        std::cout << "readInGraphRerseComplement:\n"
                  << reverseComplement << "\n";
        std::cout << "readInGraph:\n" << result << "\n";
        std::cout << "actualRead:\n" << rD->getRead(read) << "\n";
    }
    return reverseComplement == rD->getRead(read);
}

void Consensus::finishWriteConsensus(const std::vector<std::vector<read_t>>& numReadsInContig) {
    size_t size = 0;
    for (size_t i = 0; i < numThr; i++)
        size += numReadsInContig[i].size();
    std::ofstream metaData;
    metaData.open(tempDir + "metaData");
    metaData << "numReads=" << rD->getNumReads() << '\n';
    metaData << "numContigs=" << size << '\n';
    metaData << "numThr=" << numThr << '\n';
    metaData << "numReadsInContig=";
    for (size_t i = 0; i < numThr; ++i)
        for (size_t j = 0; j < numReadsInContig[i].size(); ++j)
            metaData << numReadsInContig[i][j] << ":";
    metaData << '\n';
    metaData.close();
}

ConsensusGraph *Consensus::createGraph(read_t &firstUnaddedRead) {
    // Note: firstUnaddedRead is local to the thread, and is a lower bound
    // on the actual first unadded read.
    read_t read = firstUnaddedRead;
    if (!getRead(read))
        return nullptr;
    ConsensusGraph *cG = new ConsensusGraph(aligner);
    cG->initialize(rD->getRead(read), read, 0);
    cG->calculateMainPathGreedy();
    firstUnaddedRead = read + 1;
    return cG;
}

void Consensus::initialize() {
    numReads = rD->getNumReads();
    inGraph.resize(numReads, false);
    readStatusLock.resize(numLocks);
}

bool Consensus::getRead(read_t &read) {
    if (read >= numReads)
        return false;
    while (read < numReads) {
        if (!inGraph[read]) {
            std::lock_guard<OmpMutex> lg(readStatusLock[read%numLocks]);
            // check once again inside locked region (since variables are flushed)
            if (!inGraph[read]) {
                inGraph[read] = true;
                return true;
            }
        }
        read++;
    }
    return false;
}

Consensus::Consensus() {}

