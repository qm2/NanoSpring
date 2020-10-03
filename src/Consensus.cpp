#include "Consensus.h"
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

void Consensus::generateConsensus() {
    initialize();

#pragma omp parallel
    while (hasReadsLeft()) {
        ConsensusGraph *cG = createGraph();
        ssize_t initialStartPos = cG->startPos;
        ssize_t initialEndPos = cG->endPos;
        const ssize_t len = initialEndPos - initialStartPos;
        size_t offset = len / 4;

        ssize_t curPos = cG->startPos;

        while (true) {
            std::cout << "right\n";
            addRelatedReads(cG, curPos, len);
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
            addRelatedReads(cG, curPos, len);
            curPos -= offset;
        }
    }
}

void Consensus::addRelatedReads(ConsensusGraph *cG, ssize_t curPos,
                                size_t len) {
    // Find reads likely to have overlaps
    auto stringBegin = cG->mainPath.path.begin() + curPos - cG->startPos;
    auto stringEnd =
        (ssize_t)cG->mainPath.path.size() > curPos - cG->startPos + (ssize_t)len
            ? stringBegin + len
            : cG->mainPath.path.end();
    const std::string s(stringBegin, stringEnd);
    // std::cout << curPos << "\n";
    std::vector<size_t> results;
    rF->getFilteredReads(s, results);

    // std::cout << "Found " << results.size() << " reads\n";
    // for (size_t r : results)
    //     std::cout << r << " ";
    // std::cout << '\n';

    // Try to add them one by one
    for (const size_t r : results) {
        // std::cout << "Working on read " << r << '\n';
        if (!addRead(r))
            continue;
        ssize_t relPos;
        std::string &readStr = rD->getRead(r);
        if (!rA->align(s, readStr, relPos)) {
            putReadBack(r);
            continue;
        }

        ssize_t pos = curPos + relPos;
        std::vector<Edit> editScript;
        ssize_t beginOffset, endOffset;
        if (!cG->addRead(readStr, pos, editScript, beginOffset, endOffset)) {
            putReadBack(r);
            continue;
        }
        cG->updateGraph(readStr, editScript, beginOffset, endOffset, r, pos);
        cG->calculateMainPathGreedy();
        // std::cout << "Added read " << r << " first unadded read "
        //           << firstUnaddedRead << '\n';
    }
}

void Consensus::writeConsensus() {
    size_t size = graphs.size();
    for (size_t i = 0; i < size; ++i) {
        ConsensusGraph &cG = *graphs[i];
        std::string fileName = tempFileName + std::to_string(i);
        cG.writeMainPath(fileName);
        cG.writeReads(fileName);
    }
    std::ofstream metaData;
    metaData.open(tempDir + "metaData");
    metaData << "numReads=" << rD->getNumReads() << '\n';
    metaData << "numContigs=" << size << '\n';
    metaData << "numReadsInContig=";
    for (size_t i = 0; i < size; ++i) {
        metaData << graphs[i]->getNumReads() << ":";
    }
    metaData << '\n';
    metaData.close();

    // Now we combine the files, deliminated by ".\n"
    const char *const extensions[] = {".genome", ".base", ".id", ".pos",
                                      ".type"};
    const size_t numExtensions = sizeof(extensions) / sizeof(extensions[0]);
    for (size_t i = 0; i < numExtensions; ++i) {
        std::ofstream outFile(tempDir + tempFileName + extensions[i],
                              std::ios_base::binary);
        for (size_t j = 0; j < size; j++) {
            std::string fileName =
                tempDir + tempFileName + std::to_string(j) + extensions[i];
            std::ifstream inFile(fileName, std::ios_base::binary);
            outFile << inFile.rdbuf();
            outFile << ".\n";
            inFile.close();
            boost::system::error_code ec;
            const boost::filesystem::path filePath(fileName);
            boost::filesystem::remove(filePath, ec);
        }
        outFile.close();
    }
}

ConsensusGraph *Consensus::createGraph() {
    std::lock_guard<OmpNestMutex> lg(readStatusLock);
    size_t read;
    if (!getRead(read))
        return nullptr;
    ConsensusGraph *cG = new ConsensusGraph(aligner);
    graphs.emplace_back(std::move(cG));
    cG->tempDir = tempDir;
    cG->initialize(rD->getRead(read), read, 0);
    cG->calculateMainPathGreedy();
    std::cout << "Creating graph " << graphs.size() << '\n';
    return cG;
}

void Consensus::initialize() {
    numReads = rD->getNumReads();
    inGraph.resize(numReads, false);
    firstUnaddedRead = 0;
}

bool Consensus::hasReadsLeft() {
    std::lock_guard<OmpNestMutex> lg(readStatusLock);
    return firstUnaddedRead < numReads;
}

bool Consensus::getRead(size_t &read) {
    std::lock_guard<OmpNestMutex> lg(readStatusLock);
    if (firstUnaddedRead >= numReads)
        return false;
    read = firstUnaddedRead++;
    inGraph[read] = true;
    while (firstUnaddedRead < numReads && inGraph[firstUnaddedRead])
        firstUnaddedRead++;
    return true;
}

void Consensus::putReadBack(size_t read) {
    std::lock_guard<OmpNestMutex> lg(readStatusLock);
    inGraph[read] = false;
    firstUnaddedRead = std::min(read, firstUnaddedRead);
}

bool Consensus::addRead(size_t read) {
    std::lock_guard<OmpNestMutex> lg(readStatusLock);
    if (inGraph[read])
        return false;
    inGraph[read] = true;
    while (firstUnaddedRead < numReads && inGraph[firstUnaddedRead])
        firstUnaddedRead++;
    return true;
}

Consensus::Consensus() {}

Consensus::~Consensus() {}