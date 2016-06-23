/*******************************************************************************
*                            GBS RE Site Detector                             *
*                        Copyright 2016 Kevin Murray                          *
*                        Licensed under the MPL v2.0                          *
*******************************************************************************/


#include <iostream>
#include <string>
#ifdef _OPENMP
#   include <omp.h>
#endif
#include <unordered_map>
#include <cstdio>
#include <kseq.h>
#include <zlib.h>
#include <getopt.h>

using std::unordered_map;
using std::string;
using std::cout;
using std::cerr;

KSEQ_INIT(gzFile, gzread)

int usage()
{
    cerr << "USAGE: resite -l LEN FASTQ ...\n";
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    ssize_t resite_len = 0;
    int c;
    while ((c = getopt(argc, argv, "l:")) > 0) {
        switch (c) {
            case 'l':
                resite_len = atoi(optarg);
                break;
            default:
                cerr << "bad arg " << optopt << "\n";
                usage();
                break;
        }
    }

    if (resite_len < 1 || resite_len > 100) {
        cerr << "Bad RE site length " << resite_len << "\n";
        usage();
    }
    if (optind == argc) {
        cerr << "Too few read files\n";
        usage();
    }

    unordered_map<string, size_t> mainhist;
    #pragma omp parallel for schedule(dynamic) shared(mainhist)
    for (int i = optind; i < argc; i++) {
        unordered_map<string, size_t> hist;
        const char *filename = argv[i];
        gzFile fp = gzopen(filename, "r");
        kseq_t *seq = kseq_init(fp);
        ssize_t seql = 0;
        size_t skipped = 0;
        size_t nread = 0;
        #pragma omp critical
        {
            fprintf(stderr, "Starting %s (%d of %d)\n", filename, i - optind + 1, argc - optind);
        }
        while ((seql = kseq_read(seq)) > 0) {
            nread++;
            if (seql < resite_len) {
                skipped++;
                continue;
            }
            string re(seq->seq.s, resite_len);
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
