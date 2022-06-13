// C++11 - use multiple source files.

#include <string>

#include "lest.hpp"
#include "encoding.hpp"
#include "compact.hpp"
#include "RangeMaxTree.hpp"

using namespace std;


const lest::test module[] = {

    CASE("Greedy compaction") {

        cout << "Test greedy compaction functions k=5, m=3" << endl;
        uint k = 5;
        uint m = 3;
        
        uint8_t seq[3];
        uint8_t data[2];


        SETUP( "2 kmer compaction" ) {
            cout << "\t2 kmer compaction test" << endl;
            // File with one small section prepare
            Kff_file file("compact_test.kff", "w");
            Binarizer bz(file.encoding);
            Compact comp;

            // Write needed values
            Section_GV sgv(&file);
            sgv.write_var("k", k);
            sgv.write_var("m", m);
            sgv.write_var("max", 1);
            sgv.write_var("data_size", 1);
            sgv.close();
            // Write a minimizer section
            Section_Minimizer sm(&file);
            bz.translate("AAA", 3, seq);
            sm.write_minimizer(seq);
            // First kmer
            bz.translate("TAAAG", 5, seq);
            data[0] = 1;
            sm.write_compacted_sequence(seq, 5, 1, data);
            // Second kmer
            bz.translate("AAAGC", 5, seq);
            data[0] = 2;
            sm.write_compacted_sequence(seq, 5, 0, data);
            // Second kmer
            sm.close();

            // Change from writer to reader
            file.close(false);
            file.open("r");
            sgv = Section_GV(&file);
            sgv.close();
            sm = Section_Minimizer(&file);

            SECTION( "Matrix creation tests" )
            {
                vector<vector<uint8_t *> > matrix = comp.prepare_kmer_matrix(sm);
                cout << "\t\tMatrix construction" << endl;
                EXPECT( matrix.size() == k - m + 1 );
                EXPECT( matrix[0].size() == 0u );
                
                // First kmer tests
                EXPECT( matrix[1].size() == 1u );
                EXPECT( matrix[1][0] == comp.kmer_buffer );
                bz.translate("TG", 2, seq);
                EXPECT( comp.kmer_buffer[0] == seq[0] ); // Sequence
                EXPECT( comp.kmer_buffer[1] == 1 ); // Data
                EXPECT( comp.kmer_buffer[2] == 1 ); // mini pos

                // Second kmer tests
                EXPECT( matrix[2].size() == 1u );
                bz.translate("GC", 2, seq);
                EXPECT( comp.kmer_buffer[3] == seq[0] );
                EXPECT( comp.kmer_buffer[4] == 2 );
                EXPECT( comp.kmer_buffer[5] == 0 );
            }

            sm.close();
            file.close(false);
            cout << "\t\tOK" << endl;
        }
        
        cout << endl;

    },



    CASE( "Sorted compaction" ) {

        cout << "Test sorted compaction functions k=5, m=3" << endl;
        uint k = 5;
        uint m = 3;
        
        uint8_t seq[3] = {0, 0, 0};

        SETUP( "kmer interleaved sorting" ) {
            cout << "\tkmer interleaved sorting" << endl;
            // File with one small section prepare
            // Kff_file file("compact_test.kff", "w");

            // Init the compaction object
            Compact comp;
            comp.k = k;
            comp.m = m;
            comp.data_size = 0;
            comp.mini_pos_size = 1;
            comp.bytes_compacted = 1;

            // Create the kmer matrix
            vector<vector<uint8_t *> > matrix(3, vector<uint8_t *>());

            // Add the kmers to compact
            uint8_t encoding[] = {0, 1, 3, 2};
            Binarizer bz(encoding);
            
            // AAAAA
            bz.translate("AA", k-m, seq);
            long aa_pos1 = comp.add_kmer_to_buffer(seq, nullptr, 1);
            uint8_t * aa1 = comp.kmer_buffer + aa_pos1;
            // AAAAA
            bz.translate("AA", k-m, seq);
            long aa_pos2 = comp.add_kmer_to_buffer(seq, nullptr, 2);
            uint8_t * aa2 = comp.kmer_buffer + aa_pos2;

            // GGAAA
            bz.translate("GG", k-m, seq);
            long gg_pos = comp.add_kmer_to_buffer(seq, nullptr, 2);
            uint8_t * gg = comp.kmer_buffer + gg_pos;
            // GAAC
            bz.translate("GC", k-m, seq);
            long gc_pos = comp.add_kmer_to_buffer(seq, nullptr, 1);
            uint8_t * gc = comp.kmer_buffer + gc_pos;
            // AAACT
            bz.translate("CT", k-m, seq);
            long ct_pos = comp.add_kmer_to_buffer(seq, nullptr, 0);
            uint8_t * ct = comp.kmer_buffer + ct_pos;
            // CGAAA
            bz.translate("CG", k-m, seq);
            long cg_pos = comp.add_kmer_to_buffer(seq, nullptr, 2);
            uint8_t * cg = comp.kmer_buffer + cg_pos;
            // GAAAT
            bz.translate("GT", k-m, seq);
            long gt_pos = comp.add_kmer_to_buffer(seq, nullptr, 1);
            uint8_t * gt = comp.kmer_buffer + gt_pos;
            // AAATT
            bz.translate("TT", k-m, seq);
            long tt_pos = comp.add_kmer_to_buffer(seq, nullptr, 0);
            uint8_t * tt = comp.kmer_buffer + tt_pos;
            


            SECTION( "kmer comparison" )
            {
                cout << "\t\tkmer comparison" << endl;

                // Test identity
                int cmp_ret = comp.interleaved_compare_kmers(gc, gc);
                EXPECT( cmp_ret == 0 );
                cmp_ret = comp.interleaved_compare_kmers(tt, tt);
                EXPECT( cmp_ret == 0 );
                cmp_ret = comp.interleaved_compare_kmers(gg, gg);
                EXPECT( cmp_ret == 0 );

                // Test lower than
                cmp_ret = comp.interleaved_compare_kmers(ct, tt);
                EXPECT( cmp_ret == -1 );
                cmp_ret = comp.interleaved_compare_kmers(gc, gt);
                EXPECT( cmp_ret == -1 );
                cmp_ret = comp.interleaved_compare_kmers(cg, gg);
                EXPECT( cmp_ret == -1 );

                // Test higher than
                cmp_ret = comp.interleaved_compare_kmers(tt, ct);
                EXPECT( cmp_ret == +1 );
                cmp_ret = comp.interleaved_compare_kmers(gt, gc);
                EXPECT( cmp_ret == +1 );
                cmp_ret = comp.interleaved_compare_kmers(gg, cg);
                EXPECT( cmp_ret == +1 );
            }

            matrix[0].push_back(gg);
            matrix[0].push_back(cg);
            matrix[1].push_back(gc);
            matrix[1].push_back(gt);
            matrix[2].push_back(ct);
            matrix[2].push_back(tt);

            comp.sort_matrix(matrix);
            SECTION( "Matrix sorting" )
            {
                cout << "\t\tMatrix sorting" << endl;
                // Sorting the matrix

                // First column
                EXPECT( cg == matrix[0][0] );
                EXPECT( gg == matrix[0][1] );

                // Second column
                EXPECT( gc == matrix[1][0] );
                EXPECT( gt == matrix[1][1] );

                // Third column
                EXPECT( ct == matrix[2][0] );
                EXPECT( tt == matrix[2][1] );
            }


            SECTION( "kmer pairing tests" )
            {
                cout << "\t\tkmer pairing" << endl;

                // Prepare real pairs to test
                unordered_map<uint64_t, uint64_t> real_pairs;
                //         gc   ct
                real_pairs[0] = 0;
                //         gt   tt
                real_pairs[1] = 1;

                // Perform pairing
                vector<pair<uint64_t, uint64_t> > pairs = comp.pair_kmers(matrix[1], matrix[2]);

                // Verify
                EXPECT( pairs.size() == 2u );
                for (pair<uint64_t, uint64_t> & pair : pairs) {
                    EXPECT( pair.second == real_pairs[pair.first]);
                }
            }

            matrix[0].push_back(aa1);
            matrix[1].push_back(aa2);

            SECTION( "colinear chaining test" )
            {
                cout << "\t\tBasic colinear chaining test" << endl;

                vector<pair<uint64_t, uint64_t> > pairs = comp.pair_kmers(matrix[0], matrix[1]);

                EXPECT( pairs.size() == 5u );

                // Prepare real pairs to test
                unordered_map<uint64_t, uint64_t> real_colinear;
                //             cg   gc
                real_colinear[0] = 0;
                //            gg    gt
                real_colinear[1] = 1;

                // Perform colinear chaining
                vector<pair<uint64_t, uint64_t> > co_chain = comp.colinear_chaining(pairs);

                // Verify
                EXPECT( co_chain.size() == 2u );
                // EXPECT( co_chain[0].first == cg );
                // EXPECT( co_chain[0].second == gc );
                // EXPECT( co_chain[1].first == gg );
                // EXPECT( co_chain[1].second == gt );
            }

            cout << "\t\tOK" << endl;
        }
        
        cout << endl;

    }
};

extern lest::tests & specification();

MODULE( specification(), module )