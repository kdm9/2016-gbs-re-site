/*******************************************************************************
*                            GBS RE Site Detector                             *
*                        Copyright 2016 Kevin Murray                          *
*                        Licensed under the MPL v2.0                          *
*******************************************************************************/


#include <iostream>
#include <string>
#include <omp.h>
#include <unordered_map>
#include <cstdio>
#include <kseq.h>
#include <zlib.h>

using std::unordered_map;
using std::string;
using std::cout;
using std::cerr;

KSEQ_INIT(gzFile, gzread)

#define RESITE 5 // 5 bp from start of read

int main(int argc, char *argv[])
{
    if (argc < 2) {
        cerr << "USAGE: resite FASTQ ...\n";
        return EXIT_FAILURE;
    }
    unordered_map<string, size_t> mainhist;
    #pragma omp parallel for schedule(dynamic) shared(mainhist)
    for (int i = 1; i < argc; i++) {
        unordered_map<string, size_t> hist;
        const char *filename = argv[i];
        gzFile fp = gzopen(filename, "r");
        kseq_t *seq = kseq_init(fp);
        ssize_t seql = 0;
        size_t skipped = 0;
        size_t nread = 0;
        #pragma omp critical
        {
            fprintf(stderr, "Starting %s\n", filename);
        }
        while ((seql = kseq_read(seq)) > 0) {
            nread++;
            if (seql < RESITE) {
                skipped++;
                continue;
            }
            string re(seq->seq.s, RESITE);
            hist[re] += 1;
            if (nread % 50000 == 0) {
                fprintf(stderr, ".");
                fflush(stderr);
            }
        }
        kseq_destroy(seq);
        #pragma omp critical
        {
            fprintf(stderr, "Done %s with %zu reads\n", filename, nread);
            for (auto pair: hist) {
                mainhist[pair.first] += pair.second;
            }
        }
    }
    for (auto pair: mainhist) {
        cout << pair.first << "\t" << pair.second << "\n";
    }
    return 0;
}
