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
    }
};

extern lest::tests & specification();

MODULE( specification(), module )
