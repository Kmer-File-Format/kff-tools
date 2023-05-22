// C++11 - use multiple source files.

#include "lest.hpp"
#include "skmers.hpp"


using namespace std;

const lest::test module[] = {

    CASE("Interleaved tests") {

        cout << "Test Interleaved 1 Byte" << endl;
        SETUP( "One Byte skmer" ) {
            //                   G C T A
            uint8_t skmer[] = {0b11011000};
            uint8_t inter_mem[] = {0};


            SECTION( "Minimizer position 0" )
            {
                cout << "\tMinimizer position 0" << endl;
                interleved_t inter = interleaved(skmer, inter_mem, 4, 0);

                EXPECT( inter.pref_size == 0u );
                EXPECT( inter.suf_size == 4u );
                EXPECT( inter.nucl[0] == 0b11011000 );
            }

            SECTION( "Minimizer position 1" )
            {
                cout << "\tMinimizer position 1" << endl;
                interleved_t inter = interleaved(skmer, inter_mem, 4, 1);

                EXPECT( inter.pref_size == 1u );
                EXPECT( inter.suf_size == 3u );
                EXPECT( (uint64_t)inter.nucl[0] == 0b01111000u );
            }

            SECTION( "Minimizer position 2" )
            {
                cout << "\tMinimizer position 2" << endl;
                interleved_t inter = interleaved(skmer, inter_mem, 4, 2);

                EXPECT( inter.pref_size == 2u );
                EXPECT( inter.suf_size == 2u );
                EXPECT( inter.nucl[0] == 0b10010011 );
            }

            SECTION( "Minimizer position 3" )
            {
                cout << "\tMinimizer position 3" << endl;
                interleved_t inter = interleaved(skmer, inter_mem, 4, 3);

                EXPECT( inter.pref_size == 3u );
                EXPECT( inter.suf_size == 1u );
                EXPECT( inter.nucl[0] == 0b00100111 );
            }

            SECTION( "Minimizer position 4" )
            {
                cout << "\tMinimizer position 4" << endl;
                interleved_t inter = interleaved(skmer, inter_mem, 4, 4);

                EXPECT( inter.pref_size == 4u );
                EXPECT( inter.suf_size == 0u );
                EXPECT( inter.nucl[0] == 0b00100111 );
            }

            cout << "\t\tOK" << endl;
        }

        cout << "Test Interleaved 2 Bytes" << endl;
        SETUP( "One Byte skmer" ) {
            //                   G C T A
            uint8_t skmer[] = {0b11011000, 0b11011000};
            uint8_t inter_mem[] = {0, 0};


            SECTION( "Minimizer position 3" )
            {
                cout << "\tMinimizer position 3" << endl;
                interleved_t inter = interleaved(skmer, inter_mem, 8, 3);

                EXPECT( inter.pref_size == 3u );
                EXPECT( inter.suf_size == 5u );
                EXPECT( inter.nucl[0] == 0b00101101 );
                EXPECT( inter.nucl[1] == 0b01111000 );
            }

            SECTION( "Minimizer position 4" )
            {
                cout << "\tMinimizer position 4" << endl;
                interleved_t inter = interleaved(skmer, inter_mem, 8, 4);

                EXPECT( inter.pref_size == 4u );
                EXPECT( inter.suf_size == 4u );
                EXPECT( inter.nucl[0] == 0b11000110 );
                EXPECT( inter.nucl[1] == 0b10010011 );
            }

            SECTION( "Minimizer position 5" )
            {
                cout << "\tMinimizer position 5" << endl;
                interleved_t inter = interleaved(skmer, inter_mem, 8, 5);

                EXPECT( inter.pref_size == 5u );
                EXPECT( inter.suf_size == 3u );
                EXPECT( inter.nucl[0] == 0b01111000 );
                EXPECT( inter.nucl[1] == 0b00100111 );
            }

            cout << "\t\tOK" << endl;
        }
    },

    CASE("Interleaved compare") {

        cout << "Test interleaved Compare" << endl;
        SETUP( "One Byte skmer compare" ) {
            //                   C A G T
            uint8_t skmer[] = {0b01001110};
            uint8_t skmer_test[] = {0};
            uint8_t inter_mem[] = {0};
            uint8_t inter_test_mem[] = {0};
            interleved_t inter = interleaved(skmer, inter_mem, 4, 1);
            interleved_t inter_test;

            SECTION( "Same minimizer position" )
            {
                cout << "\tSame minimizer position" << endl;

                skmer_test[0] = 0b01001110;
                inter_test = interleaved(skmer_test, inter_test_mem, 4, 1);
                EXPECT(inf_interleaved(inter, inter_test) == false);

                skmer_test[0] = 0b01111110;
                inter_test = interleaved(skmer_test, inter_test_mem, 4, 1);
                EXPECT(inf_interleaved(inter, inter_test) == true);

                skmer_test[0] = 0b10001110;
                inter_test = interleaved(skmer_test, inter_test_mem, 4, 1);
                EXPECT(inf_interleaved(inter, inter_test) == true);
                
                skmer_test[0] = 0b00001110;
                inter_test = interleaved(skmer_test, inter_test_mem, 4, 1);
                EXPECT(inf_interleaved(inter, inter_test) == false);
                
                skmer_test[0] = 0b01001100;
                inter_test = interleaved(skmer_test, inter_test_mem, 4, 1);
                EXPECT(inf_interleaved(inter, inter_test) == false);
                
                skmer_test[0] = 0b01001111;
                inter_test = interleaved(skmer_test, inter_test_mem, 4, 1);
                EXPECT(inf_interleaved(inter, inter_test) == true);
            }

            SECTION( "Different minimizer position" )
            {
                cout << "\tDifferent minimizer position" << endl;

                skmer_test[0] = 0b01001110;
                inter_test = interleaved(skmer_test, inter_test_mem, 4, 0);
                EXPECT(inf_interleaved(inter, inter_test) == true);

                skmer_test[0] = 0b01000010;
                inter_test = interleaved(skmer_test, inter_test_mem, 4, 2);
                EXPECT(inf_interleaved(inter, inter_test) == false);

                skmer_test[0] = 0b01100010;
                inter_test = interleaved(skmer_test, inter_test_mem, 4, 2);
                EXPECT(inf_interleaved(inter, inter_test) == true);
                
                skmer_test[0] = 0b00000100;
                inter_test = interleaved(skmer_test, inter_test_mem, 4, 3);
                EXPECT(inf_interleaved(inter, inter_test) == false);
                
                skmer_test[0] = 0b00000000;
                inter_test = interleaved(skmer_test, inter_test_mem, 4, 4);
                EXPECT(inf_interleaved(inter, inter_test) == false);
            }
            
            cout << "\t\tOK" << endl;
        }
    },

    CASE("Min interleaved") {

        cout << "Test minimum between interleaved" << endl;
        SETUP( "One Byte skmer minimum" ) {
            //                    TA|AA       A|AGT       AGC|T       |AGTC       ACCC|
            uint8_t skmers[] = {0b10000000, 0b00001110, 0b00110110, 0b00111001, 0b00010101};
            uint8_t inter_mems[] = {0, 0, 0, 0, 0};
            vector<interleved_t> interleaves;
            interleaves.push_back(interleaved(skmers, inter_mems, 4, 2));
            interleaves.push_back(interleaved(skmers+1, inter_mems+1, 4, 1));
            interleaves.push_back(interleaved(skmers+2, inter_mems+2, 4, 3));
            interleaves.push_back(interleaved(skmers+3, inter_mems+3, 4, 0));
            interleaves.push_back(interleaved(skmers+4, inter_mems+4, 4, 4));
            
            SECTION( "Successive minimum" )
            {
                cout << "\tSuccessive minimum tests" << endl;

                interleved_t min = min_interleaved(interleaves.begin(), interleaves.end());
                EXPECT(min.nucl[0] == 0b00000010);

                min = min_interleaved(interleaves.begin()+1, interleaves.end());
                EXPECT(min.nucl[0] == 0b00001110);

                min = min_interleaved(interleaves.begin()+2, interleaves.end());
                EXPECT(min.nucl[0] == 0b00111001);
            }
            
            cout << "\t\tOK" << endl;
        }
    }
};

extern lest::tests & specification();

MODULE( specification(), module )
