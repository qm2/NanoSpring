#ifndef F00283BB_ED79_4175_9DE1_CAEB4390DD6C
#define F00283BB_ED79_4175_9DE1_CAEB4390DD6C

#include <boost/filesystem.hpp>
#include <string>

namespace DirectoryUtils {

/**
 * @brief Combines the files filestem[0-(numFiles-1)].fileExt into the
 * single file filepath.fileExt by inplacing a ".\n" in between files, and
 * removes the old files
 *
 * @param filestem
 * @param fileExt
 * @param numFiles
 */
void combineFilesWithExt(const std::string &filestem,
                         const std::string &fileExt, size_t numFiles);

/**
 * @brief Unpacks the file filepath by splitting it up at ".\n" and stores the
 * results into outputDir/stem[0-n].ext
 *
 * @param filepath Must be a filepath having an extension
 * @param outputDir The trailing "/" is optional
 */
void unpack(const std::string &filepath, const std::string &outputDir);

/**
 * @brief Unpacks the file filepath by splitting it up at ".\n" and stores the
 * results into dir/stem[0-n].ext
 *
 * @param filepath
 */
void unpack(const std::string &filepath);

/**
 * @brief Get all extensions (.XXX) in path path and stores them in
 * extensions
 *
 * @tparam Inserter insert iterator used to store the new extensions in
 * string format.
 * @param path
 * @param extensionsIt
 */
template <typename Inserter>
void getAllExtensions(const std::string &path, Inserter extensionsIt);

/**
 * @brief Removes everything in path and creates the path directory if it does
 * not exist
 *
 * @param path
 */
void clearDir(const std::string &path);

/******************************************************************************/
/* Implementations */

template <typename Inserter>
void getAllExtensions(const std::string &path, Inserter extensionsIt) {
    boost::filesystem::directory_iterator endIt;
    for (boost::filesystem::directory_iterator it(path); it != endIt; ++it) {
        if (!boost::filesystem::is_regular_file(*it))
            continue;
        boost::filesystem::path fullPath = it->path();
        std::string filename = fullPath.filename().string();
        // We only look at files with extensions
        const auto &ext = fullPath.extension();
        if (ext.empty())
            continue;
        *(extensionsIt++) = ext.string();
    }
}

} // namespace DirectoryUtils

#endif /* F00283BB_ED79_4175_9DE1_CAEB4390DD6C */