#include "DirectoryUtils.h"
#include <boost/algorithm/string/predicate.hpp>

namespace DirectoryUtils {

void clearDir(const std::string &path) {
    boost::system::error_code ec;
    const boost::filesystem::path dirPath(path);
    boost::filesystem::remove_all(dirPath, ec);
    boost::filesystem::create_directory(dirPath, ec);
}

void combineFilesWithExt(const std::string &filestem,
                         const std::string &fileExt, const size_t numFiles) {
    std::ofstream outFile(filestem + fileExt, std::ios_base::binary);
    for (size_t i = 0; i < numFiles; i++) {
        std::string fileName = filestem + std::to_string(i) + fileExt;
        std::ifstream inFile(fileName, std::ios_base::binary);
        outFile << inFile.rdbuf();
        outFile << ".\n";
        inFile.close();
        boost::system::error_code ec;
        const boost::filesystem::path oldFilePath(fileName);
        boost::filesystem::remove(oldFilePath, ec);
    }
    outFile.close();
}

void unpack(const std::string &filepath, const std::string &outputDir) {
    std::ifstream inFile(filepath);
    boost::filesystem::path path(filepath);
    const std::string stem = path.filename().stem().string();
    const std::string ext = path.filename().extension().string();
    std::ofstream outFile;
    std::string line;
    bool need2CreateFile = true;
    size_t i = 0;
    const std::string dirAndStem =
        (boost::algorithm::ends_with(outputDir, "/") ? outputDir
                                                     : outputDir + "/") +
        stem;
    while (std::getline(inFile, line)) {
        if (need2CreateFile) {
            need2CreateFile = false;
            outFile.open(dirAndStem + std::to_string(i) + ext);
            ++i;
            outFile << line << '\n';
        } else {
            if (line == ".") {
                outFile.close();
                need2CreateFile = true;
            } else
                outFile << line << '\n';
        }
    }
}

void unpack(const std::string &filepath) {
    boost::filesystem::path path(filepath);
    unpack(filepath, path.parent_path().string());
}

} // namespace DirectoryUtils