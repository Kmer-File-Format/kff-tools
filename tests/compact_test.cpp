// C++11 - use multiple source files.

#include <string>

#include "lest.hpp"
#include "../src/encoding.hpp"
#include "../src/compact.hpp"

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

            SECTION( "Matrix creation" )
            {
                cout << "\t\tMatrix construction" << endl;
                vector<vector<long> > matrix = comp.prepare_kmer_matrix(sm);
                EXPECT( matrix.size() == k - m + 1 );
                EXPECT( matrix[0].size() == 0u );
                
                // First kmer tests
                EXPECT( matrix[1].size() == 1u );
                EXPECT( matrix[1][0] == 0l );
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


    // CASE( "Sorting tests for compaction" ) {

    //     cout << "Test sorting k=11, m=3" << endl;
    //     uint k = 5;
    //     uint m = 3;
        
    //     uint8_t seq[3] = {0, 0, 0};

    //     SETUP( "kmer interleaved sorting" ) {
    //         cout << "\tkmer interleaved sorting" << endl;
    //         // File with one small section prepare
    //         // Kff_file file("compact_test.kff", "w");

    //         // Init the compaction object
    //         Compact comp;
    //         comp.k = k;
    //         comp.m = m;
    //         comp.data_size = 0;
    //         comp.mini_pos_size = 1;
    //         comp.bytes_compacted = 1;

    //         // Create the kmer matrix
    //         vector<vector<long> > matrix(3, vector<long>());

    //         // Add the kmers to compact
    //         uint8_t encoding[] = {0, 1, 3, 2};
    //         Binarizer bz(encoding);
    //         // GGAAA
    //         bz.translate("GG", k-m, seq);
    //         long gg_pos = comp.add_kmer_to_buffer(seq, nullptr, 2);
    //         matrix[0].push_back(gg_pos);
    //         // GAAC
    //         bz.translate("GC", k-m, seq);
    //         long gc_pos = comp.add_kmer_to_buffer(seq, nullptr, 1);
    //         matrix[1].push_back(gc_pos);
    //         // AAACT
    //         bz.translate("CT", k-m, seq);
    //         long ct_pos = comp.add_kmer_to_buffer(seq, nullptr, 0);
    //         matrix[2].push_back(ct_pos);
    //         // CGAAA
    //         bz.translate("CG", k-m, seq);
    //         long cg_pos = comp.add_kmer_to_buffer(seq, nullptr, 2);
    //         matrix[0].push_back(cg_pos);
    //         // GAAAT
    //         bz.translate("GT", k-m, seq);
    //         long gt_pos = comp.add_kmer_to_buffer(seq, nullptr, 1);
    //         matrix[1].push_back(gt_pos);
    //         // AAATT
    //         bz.translate("TT", k-m, seq);
    //         long tt_pos = comp.add_kmer_to_buffer(seq, nullptr, 0);
    //         matrix[2].push_back(tt_pos);

            
    //         SECTION( "Matrix sorting" )
    //         {
    //             cout << "\t\tMatrix sorting" << endl;
    //             // Sorting the matrix
    //             comp.sort_matrix(matrix);

    //             // First column
    //             EXPECT( cg_pos == matrix[0][0] );
    //             EXPECT( gg_pos == matrix[0][1] );

    //             // Second column
    //             EXPECT( gc_pos == matrix[1][0] );
    //             EXPECT( gt_pos == matrix[1][1] );

    //             // Third column
    //             EXPECT( ct_pos == matrix[2][0] );
    //             EXPECT( tt_pos == matrix[2][1] );
    //         }

    //         cout << "\t\tOK" << endl;
    //     }
        
    //     cout << endl;

    // },


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
            vector<vector<long> > matrix(3, vector<long>());

            // Add the kmers to compact
            uint8_t encoding[] = {0, 1, 3, 2};
            Binarizer bz(encoding);
            // GGAAA
            bz.translate("GG", k-m, seq);
            long gg_pos = comp.add_kmer_to_buffer(seq, nullptr, 2);
            matrix[0].push_back(gg_pos);
            // GAAC
            bz.translate("GC", k-m, seq);
            long gc_pos = comp.add_kmer_to_buffer(seq, nullptr, 1);
            matrix[1].push_back(gc_pos);
            // AAACT
            bz.translate("CT", k-m, seq);
            long ct_pos = comp.add_kmer_to_buffer(seq, nullptr, 0);
            matrix[2].push_back(ct_pos);
            // CGAAA
            bz.translate("CG", k-m, seq);
            long cg_pos = comp.add_kmer_to_buffer(seq, nullptr, 2);
            matrix[0].push_back(cg_pos);
            // GAAAT
            bz.translate("GT", k-m, seq);
            long gt_pos = comp.add_kmer_to_buffer(seq, nullptr, 1);
            matrix[1].push_back(gt_pos);
            // AAATT
            bz.translate("TT", k-m, seq);
            long tt_pos = comp.add_kmer_to_buffer(seq, nullptr, 0);
            matrix[2].push_back(tt_pos);


            SECTION( "kmer comparison" )
            {
                cout << "\t\tkmer comparison" << endl;

                // Test identity
                int cmp_ret = comp.interleaved_compare_kmers(gc_pos, gc_pos);
                EXPECT( cmp_ret == 0 );
                cmp_ret = comp.interleaved_compare_kmers(tt_pos, tt_pos);
                EXPECT( cmp_ret == 0 );
                cmp_ret = comp.interleaved_compare_kmers(gg_pos, gg_pos);
                EXPECT( cmp_ret == 0 );

                // Test lower than
                cmp_ret = comp.interleaved_compare_kmers(ct_pos, tt_pos);
                EXPECT( cmp_ret == -1 );
                cmp_ret = comp.interleaved_compare_kmers(gc_pos, gt_pos);
                EXPECT( cmp_ret == -1 );
                cmp_ret = comp.interleaved_compare_kmers(cg_pos, gg_pos);
                EXPECT( cmp_ret == -1 );

                // Test higher than
                cmp_ret = comp.interleaved_compare_kmers(tt_pos, ct_pos);
                EXPECT( cmp_ret == +1 );
                cmp_ret = comp.interleaved_compare_kmers(gt_pos, gc_pos);
                EXPECT( cmp_ret == +1 );
                cmp_ret = comp.interleaved_compare_kmers(gg_pos, cg_pos);
                EXPECT( cmp_ret == +1 );
            }


            // SECTION( "Matrix sorting" )
            // {
            //     cout << "\t\tMatrix sorting" << endl;
            //     // Sorting the matrix
            //     comp.sort_matrix(matrix);

            //     // First column
            //     EXPECT( cg_pos == matrix[0][0] );
            //     EXPECT( gg_pos == matrix[0][1] );

            //     // Second column
            //     EXPECT( gc_pos == matrix[1][0] );
            //     EXPECT( gt_pos == matrix[1][1] );

            //     // Third column
            //     EXPECT( ct_pos == matrix[2][0] );
            //     EXPECT( tt_pos == matrix[2][1] );
            // }

            cout << "\t\tOK" << endl;
        }
        
        cout << endl;

    }
};

extern lest::tests & specification();

MODULE( specification(), module )