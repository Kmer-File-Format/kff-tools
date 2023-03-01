// C++11 - use multiple source files.

#include <string>

#include "lest.hpp"
#include "encoding.hpp"
#include "compact.hpp"
#include "RMT.hpp"

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
            vector<vector<uint8_t *> > matrix_trio(3, vector<uint8_t *>());

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

            matrix_trio[0].push_back(gg);
            matrix_trio[1].push_back(gc);
            matrix_trio[1].push_back(gt);

            comp.sort_matrix(matrix);
            comp.sort_matrix(matrix_trio);
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


            vector<vector<pair<uint64_t, uint64_t> > > pairs_by_column;
            pairs_by_column.push_back(comp.pair_kmers(matrix[0], matrix[1]));
            pairs_by_column.push_back(comp.pair_kmers(matrix[1], matrix[2]));

            vector<pair<uint64_t, uint64_t> > trio_pairs = comp.pair_kmers(matrix_trio[0], matrix_trio[1]);

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
                vector<pair<uint64_t, uint64_t> > & pairs = pairs_by_column[1];

                // Verify
                EXPECT( pairs.size() == 2u );
                for (pair<uint64_t, uint64_t> & pair : pairs) {
                    EXPECT( pair.second == real_pairs[pair.first]);
                }
            }

            vector<vector<pair<uint64_t, uint64_t> > > colinear_chainings;
            for (const auto & pair : pairs_by_column)
                colinear_chainings.push_back(comp.colinear_chaining(pair));

            SECTION( "colinear chaining test 1" )
            {
                cout << "\t\tBasic colinear chaining test 1" << endl;
                matrix[0].push_back(aa1);
                matrix[1].push_back(aa2);

                vector<pair<uint64_t, uint64_t> > pairs = comp.pair_kmers(matrix[0], matrix[1]);

                EXPECT( pairs.size() == 5u );

                // Prepare real pairs to test
                unordered_map<uint64_t, uint64_t> real_colinear;
                //           cg    gc
                real_colinear[0] = 0;
                //           gg    gt
                real_colinear[1] = 1;
                //           aa    aa
                real_colinear[2] = 2;

                // Perform colinear chaining
                vector<pair<uint64_t, uint64_t> > co_chain = comp.colinear_chaining(pairs);

                // Verify
                EXPECT( co_chain.size() == 3u );
                for (pair<uint64_t, uint64_t> & p : co_chain)
                    EXPECT( real_colinear[p.first] == p.second );

                cout << "\t\tSecond colinear chaining test" << endl;

                pairs = {{0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 0}};

                // Perform colinear chaining
                co_chain = comp.colinear_chaining(pairs);

                // Verify
                EXPECT( co_chain.size() == 2u );
                EXPECT( co_chain[0].first == 0u); EXPECT(co_chain[0].second == 1u);
                EXPECT( co_chain[1].first == 1u); EXPECT(co_chain[1].second == 2u);

                /// Other case

                pairs = {{0, 1}, {0, 2}, {1, 1}, {2, 0}, {2, 1}};

                // Perform colinear chaining
                co_chain = comp.colinear_chaining(pairs);

                // Verify
                EXPECT( co_chain.size() == 1u);
                EXPECT( co_chain[0].first == 0u); EXPECT(co_chain[0].second == 2u);

                /// One more case

                pairs = {{0, 1}, {0, 2}, {1, 0}, {1, 1}, {2, 1}};

                // Perform colinear chaining
                co_chain = comp.colinear_chaining(pairs);

                // Verify
                EXPECT( co_chain.size() == 2u);
                EXPECT( co_chain[0].first == 1u); EXPECT(co_chain[0].second == 0u);
                EXPECT( co_chain[1].first == 2u); EXPECT(co_chain[1].second == 1u);

                /// One more case

                pairs = {{0, 0}, {0, 1}};

                // Perform colinear chaining
                co_chain = comp.colinear_chaining(pairs);

                // Verify
                EXPECT( co_chain.size() == 1u);
                EXPECT( co_chain[0].first == 0u); EXPECT(co_chain[0].second == 1u);

                /// One more case

                pairs = {{0, 0}, {1, 0}};

                // Perform colinear chaining
                co_chain = comp.colinear_chaining(pairs);

                // Verify

                EXPECT( co_chain.size() == 1u);
                EXPECT( co_chain[0].first == 1u); EXPECT(co_chain[0].second == 0u);

                /// One more case

                pairs = {{0, 1}, {1, 0}};

                // Perform colinear chaining
                co_chain = comp.colinear_chaining(pairs);

                // Verify
                EXPECT( co_chain.size() == 1u);
                EXPECT( co_chain[0].first == 0u); EXPECT(co_chain[0].second == 1u);

                /// One more case

                pairs = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};

                // Perform colinear chaining
                co_chain = comp.colinear_chaining(pairs);

                // Verify
                EXPECT( co_chain.size() == 2u);
                EXPECT( co_chain[0].first == 0u); EXPECT(co_chain[0].second == 0u);
                EXPECT( co_chain[1].first == 1u); EXPECT(co_chain[1].second == 1u);


            }


            SECTION( "colinear chaining test 2" )
            {
                cout << "\t\tBasic colinear chaining test 2" << endl;

                EXPECT( trio_pairs.size() == 2u );

                // Perform colinear chaining
                vector<pair<uint64_t, uint64_t> > trio_chain = comp.colinear_chaining(trio_pairs);

                // Verify
                EXPECT( trio_chain.size() == 1u );
                EXPECT( 1u == trio_chain[0].second ); // 0 needed ?
            }

            SECTION(" colinear chaining test 3") {
                cout << "\t\tBasic colinear chaining test 3" << endl;

                // Making pairs
                vector<pair<uint64_t, uint64_t> > pairs = {{0, 0}, {0, 1}, {1, 1}, {1, 2}, {1, 3}};

                // Perform colinear chaining
                vector<pair<uint64_t, uint64_t> > co_chain = comp.colinear_chaining(pairs);

                // Verify
                EXPECT( co_chain.size() == 2u);
                EXPECT( co_chain[0].first == 0u); EXPECT(co_chain[0].second == 0u);
                EXPECT( co_chain[1].first == 1u); EXPECT(co_chain[1].second == 3u);

                /// Reverse exemple

                // Making pairs
                pairs = {{0, 0}, {1, 1}, {2, 1}, {3, 1}};

                // Perform colinear chaining
                co_chain = comp.colinear_chaining(pairs);

                // Verify

                EXPECT( co_chain.size() == 2u);
                EXPECT( co_chain[0].first == 0u); EXPECT(co_chain[0].second == 0u);
                EXPECT( co_chain[1].first == 3u); EXPECT(co_chain[1].second == 1u);
            }


            SECTION( "polish super-kmers" )
            {
                cout << "\t\tPolish super-kmers test" << endl;

                vector<vector<uint8_t *> > true_skmers;
                true_skmers.push_back({cg, gc, ct});
                true_skmers.push_back({gg, gt, tt});

                vector<vector<uint8_t *> > skmers = comp.polish_sort(matrix, colinear_chainings);

                // Verify
                EXPECT( skmers.size() == true_skmers.size() );
                for (size_t sk_idx=0 ; sk_idx<true_skmers.size() ; sk_idx++) {
                    vector<uint8_t *> & skmer = skmers[sk_idx];
                    vector<uint8_t *> & true_sk = true_skmers[sk_idx];

                    EXPECT(skmer.size() == true_sk.size());

                    for (size_t kmer_idx=0 ; kmer_idx<true_sk.size() ; kmer_idx++)
                        EXPECT( skmer[kmer_idx] == true_sk[kmer_idx] );
                }
            }

            cout << "\t\tOK" << endl;
        }
        
        cout << endl;

    }
};

extern lest::tests & specification();

MODULE( specification(), module )