#ifndef BB85B7EA_F0A5_4EEB_BC52_352F8C5D617B
#define BB85B7EA_F0A5_4EEB_BC52_352F8C5D617B

#include "Edits.h"
#include "LocalMyers.h"
#include "MyersAligner.h"

/***
 * This class basically runs the local myers algorithm from left the right to
 * obtain a good alignment at the right most end, and then uses that alignment
 * to run the local myers algorithm from right to left. The old version that
 * doesn't consider substitutions.
 */
class LocalMyersRollBackOld : public LocalMyers {
public:
    /***
     * The maximum edit distance to search for among the overlapping portions
     */
    const size_t maxEditDis;

    const double errorRate = 0.2;

    virtual bool align(const std::string &s1, const std::string &s2,
                       std::vector<Edit> &editScript, size_t &editDis) override;

    virtual bool align(const std::string &s1, const std::string &s2,
                       const ssize_t offsetGuess, ssize_t &beginOffset,
                       ssize_t &endOffset, std::vector<Edit> &editScript,
                       size_t &editDis) override;

    LocalMyersRollBackOld(const size_t lenA, const size_t lenB,
                          const size_t maxEditDis);

    static bool localAlign(const char *&Abegin, const char *const Aend,
                           const char *&Bbegin, const char *const Bend,
                           const size_t max, std::vector<Edit> &editScript,
                           size_t &editDis, bool forward);

private:
    bool alignReverse(const std::string &s1, const std::string &s2,
                      const ssize_t offsetGuess, ssize_t &beginOffset,
                      ssize_t &endOffset, std::vector<Edit> &editScript,
                      size_t &editDis);
};

#endif /* BB85B7EA_F0A5_4EEB_BC52_352F8C5D617B */
