// C++11 - use multiple source files.

#include <string>

#include "lest.hpp"
#include "../src/encoding.hpp"
#include "../src/sequences.hpp"

using namespace std;


const lest::test module[] = {

    CASE("Testing sequence minimizer search") {

        cout << "Test minimizer search\n\tk=11, m=3, single side minimizer search" << endl;
        SETUP( "k=11, m=3, single side minimizer search" ) {
            uint8_t encoding[] = {0, 1, 3, 2};
            uint k = 5;
            uint m = 3;
            bool single_side = true;
            uint max_seq_size = 3 * k;
            uint max_nb_mmers = max_seq_size - m + 1;
            uint max_nb_kmers = max_seq_size - k + 1;

            MinimizerSearcher ms(k, m, encoding, max_seq_size, single_side);
            Binarizer bz(encoding);
            
            string seq("CATTGCA");
            uint8_t bin[2];
            std::string nucleotides(seq);
            bz.translate(nucleotides, seq.length(), bin);

            uint seq_nb_mmers = seq.length() - m + 1;
            uint seq_nb_kmers = seq.length() - k + 1;

            SECTION( "Constructor" )
            {
                cout << "\t\tConstructor" << endl;
                EXPECT( ms.mini_buffer.size() == 2 * max_nb_mmers );
                EXPECT( ms.mini_pos.size() == max_nb_kmers );
            } 

            SECTION( "Candidates generation" )
            {
                cout << "\t\tCandidates generation" << endl;
                                            //     CAT       ATT       TTG       TGC       GCA
                uint64_t awaited_candidates[] = {0b010010, 0b001010, 0b101011, 0b101101, 0b110100};

                ms.compute_candidates(bin, seq.length());
                for (uint idx=0 ; idx<seq_nb_mmers ; idx++) {
                    EXPECT( ms.mini_buffer[idx] == awaited_candidates[idx]);
                }
            }

            SECTION( "Minimizer computation" )
            {
                cout << "\t\tMinimizers generation" << endl;
                                     //   ATT  ATT  TTG   
                int64_t awaited_pos[] = {1  , 1  , 2};

                ms.compute_candidates(bin, seq.length());
                ms.compute_minimizers(seq_nb_kmers);
                for (uint idx=0 ; idx<seq_nb_kmers ; idx++) {
                    EXPECT( ms.mini_pos[idx] == awaited_pos[idx]);
                }
            }

            SECTION( "Skmer borders computation" )
            {
                cout << "\t\tSkmers delimitations" << endl;
                                     //   ATT  ATT  TTG   
                pair<uint64_t, uint64_t> awaited_skmers[] = {make_pair(0, k), make_pair(2, k+1)};

                ms.compute_candidates(bin, seq.length());
                ms.compute_minimizers(seq_nb_kmers);
                ms.compute_skmers(seq_nb_kmers);

                EXPECT( ms.skmers.size() == 2u );

                for (uint idx=0 ; idx<ms.skmers.size() ; idx++) {
                    EXPECT( ms.skmers[idx] == awaited_skmers[idx]);
                }
            }

            SECTION( "Skmer get" )
            {
                cout << "\t\tSkmers list" << endl;

                skmer awaited_skmers[] = {{0, k, 1, 0b001010}, {2, k+1, 2, 0b101011}};
                vector<skmer> skmers = ms.get_skmers(bin, seq.length());

                EXPECT( skmers.size() == 2u );

                for (uint idx=0 ; idx<skmers.size() ; idx++) {
                    EXPECT( skmers[idx].start_position == awaited_skmers[idx].start_position);
                    EXPECT( skmers[idx].stop_position == awaited_skmers[idx].stop_position);
                    EXPECT( skmers[idx].minimizer_position == awaited_skmers[idx].minimizer_position);
                    EXPECT( skmers[idx].minimizer == awaited_skmers[idx].minimizer);
                }
            }

            cout << "\t\tOK" << endl;
        }
        cout << endl;

        cout << "\tk=11, m=3, double side minimizer search" << endl;
        SETUP( "k=11, m=3, double side minimizer search" ) {
            uint8_t encoding[] = {0, 1, 3, 2};
            uint k = 5;
            uint m = 3;
            bool single_side = false;
            uint max_seq_size = 3 * k;
            uint max_nb_mmers = max_seq_size - m + 1;
            uint max_nb_kmers = max_seq_size - k + 1;

            MinimizerSearcher ms(k, m, encoding, max_seq_size, single_side);
            Binarizer bz(encoding);
            
            string seq("CATTGCA");
            uint8_t bin[2];
            std::string nucleotides(seq);
            bz.translate(nucleotides, seq.length(), bin);
            
            uint seq_nb_mmers = seq.length() - m + 1;
            uint seq_nb_kmers = seq.length() - k + 1;

            SECTION( "Constructor" )
            {
                cout << "\t\tConstructor" << endl;
                EXPECT( ms.mini_buffer.size() == 2 * max_nb_mmers );
                EXPECT( ms.mini_pos.size() == max_nb_kmers );
            } 

            SECTION( "Candidates generation" )
            {
                cout << "\t\tCandidates generation [reverse]" << endl;
                                            //     ATG       AAT       CAA       GCA       TGC
                uint64_t awaited_candidates[] = {0b001011, 0b000010, 0b010000, 0b110100, 0b101101};
                
                ms.compute_candidates(bin, seq.length());
                for (uint idx=0 ; idx<seq_nb_mmers ; idx++) {
                    EXPECT( ms.mini_buffer[max_nb_mmers + idx] == awaited_candidates[idx]);
                }
            }

            SECTION( "Minimizer computation" )
            {
                cout << "\t\tMinimizers generation [reverse]" << endl;
                                     //  AAT AAT TTG   
                int64_t awaited_pos[] = {-2, -2, -3};

                ms.compute_candidates(bin, seq.length());
                ms.compute_minimizers(seq_nb_kmers);
                for (uint idx=0 ; idx<seq_nb_kmers ; idx++) {
                    EXPECT( ms.mini_pos[idx] == awaited_pos[idx]);
                }
            }

            SECTION( "Skmer get" )
            {
                cout << "\t\tSkmers list" << endl;

                skmer awaited_skmers[] = {{0, k, -2, 0b000010}, {2, k+1, -3, 0b010000}};
                vector<skmer> skmers = ms.get_skmers(bin, seq.length());

                EXPECT( skmers.size() == 2u );

                for (uint idx=0 ; idx<skmers.size() ; idx++) {
                    EXPECT( skmers[idx].start_position == awaited_skmers[idx].start_position);
                    EXPECT( skmers[idx].stop_position == awaited_skmers[idx].stop_position);
                    EXPECT( skmers[idx].minimizer_position == awaited_skmers[idx].minimizer_position);
                    EXPECT( skmers[idx].minimizer == awaited_skmers[idx].minimizer);
                }
            }

            cout << "\t\tOK" << endl;
        }

    }
};

extern lest::tests & specification();

MODULE( specification(), module )
