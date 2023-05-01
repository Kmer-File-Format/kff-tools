// C++11 - use multiple source files.
#include "lest.hpp"
#include "colinear.hpp"


using namespace std;


const lest::test module[] = {

    CASE("Test colinear chaining on kmers") {
        cout << "kmer colinear chaining tests" << endl;
        
        SETUP("Basic construction tests") {

            SECTION( "Construction basic" )
            {
                vector<PairInt> pairs({PairInt(0, 0), PairInt(2, 2), PairInt(1, 1), PairInt(3, 3)});
                Colinear col(pairs);
                cout << "\tConstruction test" << endl;
                vector<PairInt> expected_keys{
                    PairInt(0, 0), PairInt(1, 1), PairInt(1, 1), PairInt(3, 3),
                    PairInt(2, 2), PairInt(3, 3), PairInt(3, 3)
                };

                EXPECT( col.nc_tree.size() == expected_keys.size() );

                for (uint i=0 ; i<expected_keys.size() ; i++) {
                    EXPECT( col.nc_tree[i].first == expected_keys[i] );
                    EXPECT( col.nc_tree[i].second == 0u );
                }
                cout << "\t\tOK" << endl;
            }

            SECTION( "Construction cross" )
            {
                cout << "\tConstruction cross test" << endl;
                vector<PairInt> pairs({PairInt(0, 0), PairInt(2, 1), PairInt(1, 2), PairInt(3, 3)});
                Colinear col(pairs);
                vector<PairInt> expected_keys{
                    PairInt(0, 0), PairInt(2, 1), PairInt(2, 1), PairInt(3, 3),
                    PairInt(1, 2), PairInt(3, 3), PairInt(3, 3)
                };

                EXPECT( col.nc_tree.size() == expected_keys.size() );

                for (uint i=0 ; i<expected_keys.size() ; i++) {
                    EXPECT( col.nc_tree[i].first == expected_keys[i] );
                    EXPECT( col.nc_tree[i].second == 0u );
                }
                cout << "\t\tOK" << endl;
            }
        }

        SETUP("4 pairs basic tests") {
            cout << "\t4 pairs tests" << endl;

            SECTION( "No collision" ) {
                cout << "\t\tNo collision" << endl;
                vector<PairInt> pairs(
                    {PairInt(0, 0), PairInt(2, 2), PairInt(1, 1), PairInt(3, 3)}
                );
                vector<PairInt> expected_chain(
                    {PairInt(0, 0), PairInt(1, 1), PairInt(2, 2), PairInt(3, 3)}
                );

                Colinear col(pairs);
                vector<PairInt> chain = col.longest_chain();
                
                EXPECT(chain.size() == expected_chain.size());

                for (uint64_t i(0) ; i<chain.size() ; i++)
                    EXPECT(chain[i] == expected_chain[i]);

                cout << "\t\tOK" << endl;
            }

            SECTION( "Middle collision" ) {
                cout << "\t\tMiddle collision" << endl;
                vector<PairInt> pairs(
                    {PairInt(0, 0), PairInt(2, 1), PairInt(1, 2), PairInt(3, 3)}
                );
                vector<PairInt> expected_chain(
                    {PairInt(0, 0), PairInt(1, 2), PairInt(3, 3)}
                );

                Colinear col(pairs);
                vector<PairInt> chain = col.longest_chain();
                
                EXPECT(chain.size() == expected_chain.size());

                for (uint64_t i(0) ; i<chain.size() ; i++)
                    EXPECT(chain[i] == expected_chain[i]);

                cout << "\t\tOK" << endl;
            }

            SECTION( "External collision" ) {
                cout << "\t\tMiddle collision" << endl;
                vector<PairInt> pairs(
                    {PairInt(0, 3), PairInt(2, 2), PairInt(1, 1), PairInt(3, 0)}
                );
                vector<PairInt> expected_chain(
                    {PairInt(1, 1), PairInt(2, 2)}
                );

                Colinear col(pairs);
                vector<PairInt> chain = col.longest_chain();
                
                EXPECT(chain.size() == expected_chain.size());

                for (uint64_t i(0) ; i<chain.size() ; i++)
                    EXPECT(chain[i] == expected_chain[i]);

                cout << "\t\tOK" << endl;
            }
        }
    }
};

extern lest::tests & specification();

MODULE( specification(), module )
