// C++11 - use multiple source files.

#include <string>

#include "lest.hpp"
#include "../src/encoding.hpp"
#include "../src/compact.hpp"

using namespace std;


const lest::test module[] = {

    CASE("Gredy compaction") {

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

    }
};

extern lest::tests & specification();

MODULE( specification(), module )