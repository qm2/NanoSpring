#include "Consensus.h"
#include "DirectoryUtils.h"
#include "bsc_helper.h"
#include "minimap.h"
#include <arpa/inet.h>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <csignal>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <set>
#include <ctime>
#include <chrono>
#include <malloc.h>

void Consensus::generateAndWriteConsensus() {
    initialize();

    std::vector<std::vector<read_t>> numReadsInContig(numThr);
    std::vector<std::vector<read_t>> loneReads(numThr);
    ConsensusGraph *cG = nullptr;
    std::vector<CountStats> count_stats(numThr);

#pragma omp parallel private(cG)
    {
        auto tid = omp_get_thread_num();   
        std::ofstream logfile;
#ifdef LOG
        logfile.open("logfile"+std::to_string(tid), std::ofstream::out);
#endif
        std::string filePrefix = tempDir + tempFileName + ".tid." + std::to_string(tid);
        ConsensusGraphWriter cgw(filePrefix);
        read_t firstUnaddedRead = 0;
        int contigId = 0;
        // guarantee that all reads < firstUnaddedRead have been picked
        while ((cG = createGraph(firstUnaddedRead))) {
#ifdef LOG
            logfile<< "Contig: " << contigId << ", First read number "<< cG->firstReadId <<", read length: " << cG->endPos - cG->startPos << "\n";
            {
            auto end = std::chrono::system_clock::now();
            std::time_t end_time = std::chrono::system_clock::to_time_t(end);
            logfile <<"Time: "<<std::ctime(&end_time);
            }
#endif 
            ssize_t initialStartPos = cG->startPos; // simply 0
            ssize_t initialEndPos = cG->endPos; // 
            const ssize_t len = initialEndPos - initialStartPos;
            // TODO: Experiment with this
            size_t offset = rD->avgReadLen / 4;

            ssize_t curPos = cG->startPos;
            // flag for too many edges in the graph
            bool edgesTooMany = false;
            while (len >= 32) { // for len below 32, minhash is meaningless
#ifdef LOG
                std::cout << "right\n";
#endif
                addRelatedReads(cG, curPos, len, count_stats[tid], logfile, contigId);
#ifdef LOG
                cG->printStatus();
#endif
                curPos += offset;
                // std::cout << "curPos " << curPos << " len " << len << " endPos "
                //          << cG->endPos << '\n';

                if (curPos + len > cG->endPos){  
                    break;
                }else if (cG->getNumEdges()>=edge_threshold){
                    edgesTooMany = true;
                    break;
                }
            }

            curPos = initialStartPos - offset;
            while ((len >= 32) && (!edgesTooMany)) { // for len below 32, minhash is meaningless
#ifdef LOG
                std::cout << "left\n";
#endif
                if (curPos < cG->startPos){
                    break;
                } else if (cG->getNumEdges()>=edge_threshold){  
                    edgesTooMany = true;
                    break;
                }
#ifdef LOG
                cG->printStatus();
#endif
                addRelatedReads(cG, curPos, len, count_stats[tid], logfile, contigId);
                curPos -= offset;
            }

            // if the graph has just one read, write to lone
            if (cG->getNumReads() == 0) {
                cG->writeReadLone(cgw);
                loneReads[tid].push_back(cG->firstReadId);
                numReadsInContig[tid].push_back(1);
            } else {
                cG->writeMainPath(cgw);
                cG->writeReads(cgw);
                numReadsInContig[tid].push_back(cG->getNumReads());
            }
            // if the graph is large, run malloc_trim so that memory released to system
            // without this the memory keeps increasing to very high
            bool run_malloc_trim = false;
            if (cG->getNumEdges() > 1000000)
                run_malloc_trim = true;

            delete cG;
            if (run_malloc_trim)
                malloc_trim(0);

            contigId++;
#ifdef LOG
            {
            auto end = std::chrono::system_clock::now();
            std::time_t end_time = std::chrono::system_clock::to_time_t(end);
            logfile <<"Time: "<<std::ctime(&end_time);
            }
#endif 
        }
        // finally, write ids of lone reads to end of id file
        cG->writeIdsLone(cgw, loneReads[tid]);
#ifdef LOG
        {
        auto end = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        logfile<<"Thread "<<tid<<" ends at time: "<<std::ctime(&end_time);
        }
        logfile.close();
#endif
    } // pragma omp parallel

    // now perform last step, combining files from threads and writing metadata
    finishWriteConsensus(numReadsInContig);

    // print stats about lone reads
    read_t totalNumLoneReads = 0;
    std::cout << "LoneReads";
    for (auto &v: loneReads) {
        totalNumLoneReads += v.size();
#ifdef LOG
        for (auto &r: v)
           std::cout << std::dec << ":" << r; 
#endif
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
#ifdef LOG
            logfile<<"Contig: " << contigId << ", Found "<< results.size() << " MinHash results\n";
#endif 
        cs.countMinHash += results.size();

        // Try to add them one by one
        for (const auto r : results) {
        	check if we exceed the edge limit in the graph
        	if (cG->getNumEdges()>=edge_threshold){
            	return;
        	}
            if (inGraph[r])
                continue;

            cs.countMinHashNotInGraph++;

            std::string readStr, readStr1;
            rD->getRead(r,readStr1);
            if (readStr1.size() < 32){ // for len below 32, minhash is meaningless
                continue;
            }
#ifdef LOG
            logfile<<"Contig: " << contigId << ", Read passed MinHash "<<r<<", read length: " << readStr1.length()<< "\n";
#endif 
            if (reverseComplement)
                ReadData::toReverseComplement(
                    readStr1.begin(), readStr1.end(),
                    std::inserter(readStr, readStr.end()));
            else
                readStr = readStr1;

            bool use_sort_merge = false;
            ssize_t pos_sort_merge;
            if (use_sort_merge) {
                ssize_t relPos;            
                if (!rA->align(originalString, readStr, relPos)) {
#ifdef LOG
                    logfile<< "Contig: " << contigId << ", Read failed Sort-Merge "<<r<<"\n";
#endif 
                    continue;
                }
#ifdef LOG
                logfile<<"Contig: " << contigId << ", Read passed Sort-Merge "<<r<<"\n";
#endif
                cs.countMergeSort++;
                pos_sort_merge = curPos + relPos;
            }

            std::vector<Edit> editScript;
            ssize_t beginOffset, endOffset;
            ssize_t pos;
            bool alignStatus = cG->alignRead(readStr, editScript, pos, beginOffset,
                                endOffset, m_k, m_w, max_chain_iter);
            
            if (!alignStatus) {
#ifdef LOG
                logfile<< "Contig: " << contigId << ", Read failed aligner "<<r<<"\n";
#endif
                continue;
            } else {
                if (!readStatusLock[r%numLocks].try_lock()) {
                    // we only try_lock here since missing a read
                    // is not a major issue and lock contention should be a rare event anyway.
                    // Note that if some other thread has locked the read, they are guaranteed
                    // to pick it up.
                    continue;
                } else {
                    // check again that read is available (variables flushed after lock is set)
                    if (inGraph[r]) {
                        // read already taken, continue with next read
                        readStatusLock[r%numLocks].unlock();
                        continue;
                    }

                    // read added to graph
#ifdef LOG
                    logfile<< "Contig: " << contigId << ", Read passed aligner "<<r<<"\n";
#endif
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
            // if the graph is not initialized, first initialize
            if (cG->getNumReads() == 0) {
                std::string mainPathString(cG->mainPath.path);
                cG->mainPath.path.clear();
                cG->initialize(mainPathString, cG->firstReadId, 0);
                cG->calculateMainPathGreedy();
            }
            cG->updateGraph(readStr, editScript, beginOffset, endOffset, r, pos,
                            reverseComplement);
#ifdef CHECKS
            assert(checkRead(cG, r));
            assert(cG->checkNoCycle());
#endif
            cG->calculateMainPathGreedy();
#ifdef CHECKS
            // std::cout << "Added read " << r << " first unadded read "
            //           << firstUnaddedRead << std::endl;
            assert(checkRead(cG, r));
            assert(cG->checkNoCycle());
#endif
        }
    }
}

bool Consensus::checkRead(ConsensusGraph *cG, read_t read) {
    std::string result;
    std::string readStr;
    bool temp = cG->getRead(read, std::inserter(result, result.end()));
    assert(temp);
    if (!temp)
        return false;
    rD->getRead(read, readStr);
    if (!cG->readsInGraph.at(read).reverseComplement) {
        if (result != readStr) {
            std::cout << "readInGraph:\n" << result << "\n";
            std::cout << "actualRead:\n" << readStr << "\n";
        }
        return result == readStr;
    }
    std::string reverseComplement;
    ReadData::toReverseComplement(
        result.begin(), result.end(),
        std::inserter(reverseComplement, reverseComplement.end()));
    if (reverseComplement != readStr) {
        std::cout << "readInGraphRerseComplement:\n"
                  << reverseComplement << "\n";
        std::cout << "readInGraph:\n" << result << "\n";
        std::cout << "actualRead:\n" << readStr << "\n";
    }
    return reverseComplement == readStr;
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
    ConsensusGraph *cG = new ConsensusGraph;
    // we don't actually fully initialize the graph
    // since that would be wasted effort if this turns out to be a lone graph
    rD->getRead(read, cG->mainPath.path);
    cG->startPos = 0;
    cG->endPos = cG->mainPath.path.size();
    cG->firstReadId = read;
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
            if (!readStatusLock[read%numLocks].try_lock()) {
                // couldn't obtain lock, that's fine, some other thread
                // will take the read
                read++;
            } else {
                // check once again inside locked region (since variables are flushed)
                if (!inGraph[read]) {
                    inGraph[read] = true;
                    readStatusLock[read%numLocks].unlock();
                    return true;
                }
                readStatusLock[read%numLocks].unlock();
                read++;
            }
        } else {
            read++;
        }
    }
    return false;
}

Consensus::Consensus() {}

