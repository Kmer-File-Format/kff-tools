// C++11 - use multiple source files.

#include <string>

#include "lest.hpp"
#include "../src/encoding.hpp"



const lest::test module[] = {

    CASE("A binarizer can translate strings of different sizes") {

        SETUP( "Binarizer on basic encoding" ) {
            uint8_t encoding[] = {0, 1, 3, 2};
            Binarizer bz(encoding);

            SECTION( "Test 2 Bytes bin" )
            {
                std::string nucleotides("CATTAGCA");
                uint8_t bin[2];
                bz.translate(nucleotides, bin);

                EXPECT( (uint)bin[0] == (uint)0b01001010 );
                EXPECT( (uint)bin[1] == (uint)0b00110100 );


                nucleotides = "TCGTACG";
                bz.translate(nucleotides, bin);
                
                EXPECT( (uint)bin[0] == (uint)0b100111 );
                EXPECT( (uint)bin[1] == (uint)0b10000111 );


                nucleotides = "CGTACG";
                bz.translate(nucleotides, bin);
                
                EXPECT( (uint)bin[0] == (uint)0b0111 );
                EXPECT( (uint)bin[1] == (uint)0b10000111 );


                nucleotides = "GTACG";
                bz.translate(nucleotides, bin);
                
                EXPECT( (uint)bin[0] == (uint)0b11 );
                EXPECT( (uint)bin[1] == (uint)0b10000111 );
            }
        }
    }
};

extern lest::tests & specification();

MODULE( specification(), module )
