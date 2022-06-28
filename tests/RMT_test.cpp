// C++11 - use multiple source files.

#include "lest.hpp"
#include "RMT.hpp"


using namespace std;

const lest::test module[] = {

    CASE("Test a 0-Range Max Tree implementation") {
        cout << "Test 0-Range Max Tree" << endl;
        
        SETUP("Init") {
            vector<uint64_t> keys{0, 1, 2, 3, 4, 5, 6, 7};
            RangeMaxTree<uint64_t> rmt(keys);

            SECTION( "Test 2 Bytes bin" )
            {
                cout << "  Construction" << endl;
                EXPECT( rmt.tree.size() == 15u );

                vector<uint64_t> expected_keys{0, 1, 1, 3, 2, 3, 3, 7, 4, 5, 5, 7, 6, 7, 7};
                for (uint i=0 ; i<15 ; i++) {
                    EXPECT( rmt.tree[i].first == expected_keys[i] );
                    EXPECT( rmt.tree[i].second == 0u );
                }
                cout << "  OK" << endl;
            }

            SECTION( "Update" ) {
                cout << "  Update" << endl;

                rmt.update(3, 2);
                vector<uint64_t> expected_values1{0, 0, 0, 2, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0};
                for (uint idx=0 ; idx<rmt.tree.size() ; idx++)
                    EXPECT(rmt.tree[idx].second == expected_values1[idx]);

                rmt.update(5, 3);
                vector<uint64_t> expected_values2{0, 0, 0, 2, 0, 2, 2, 3, 0, 3, 3, 3, 0, 0, 0};
                for (uint idx=0 ; idx<rmt.tree.size() ; idx++)
                    EXPECT(rmt.tree[idx].second == expected_values2[idx]);
                cout << "  OK" << endl;
            }

            SECTION( "0-Range Max" ) {
                cout << "  0-Range Max" << endl;
                
                rmt.update(3, 2);
                rmt.update(5, 3);

                EXPECT( rmt.zero_range(1) == 0u );
                EXPECT( rmt.zero_range(3) == 2u );
                EXPECT( rmt.zero_range(4) == 2u );
                EXPECT( rmt.zero_range(6) == 3u );

                EXPECT( rmt.first_max_key(0) == 0u );
                EXPECT( rmt.first_max_key(2) == 3u );
                EXPECT( rmt.first_max_key(3) == 5u );

                cout << "  OK" << endl;
            }

            SECTION( "Range Max" ) {
                cout << "  Range Max" << endl;
                
                rmt.update(1, 7);
                rmt.update(2, 2);
                rmt.update(3, 9);
                rmt.update(5, 10);

                EXPECT( rmt.range(1, 2) == 7u );
                EXPECT( rmt.range(0, 2) == 7u );
                EXPECT( rmt.range(2, 2) == 2u );
                EXPECT( rmt.range(2, 4) == 9u );
                EXPECT( rmt.range(6, 7) == 0u );

                cout << "  OK" << endl;
            }

            SECTION( "Bounded first max" ) {
                cout << "  Bounded first max" << endl;
                
                rmt.update(1, 7);
                rmt.update(2, 2);
                rmt.update(3, 9);
                rmt.update(5, 10);
                rmt.update(6, 7);

                EXPECT( rmt.bounded_first_max_key(9, 1) == 3u );
                EXPECT( rmt.bounded_first_max_key(10, 1) == 5u );
                EXPECT( rmt.bounded_first_max_key(7, 1) == 1u );
                EXPECT( rmt.bounded_first_max_key(7, 2) == 3u );

                cout << "  OK" << endl;
            }
        }

    }
};

extern lest::tests & specification();

MODULE( specification(), module )
