// C++11 - use multiple source files.

#include "lest.hpp"
#include "ZRMT.hpp"


using namespace std;

const lest::test module[] = {

    CASE("Test a 0-Range Max Tree implementation") {
        cout << "Test 0-Range Max Tree" << endl;
        
        SETUP("Init") {
            vector<uint64_t> keys{0, 1, 2, 3, 4, 5, 6, 7};
            ZeroRangeMaxTree<uint64_t> zrmt(keys);

            SECTION( "Test 2 Bytes bin" )
            {
                cout << "  Construction" << endl;
                EXPECT( zrmt.tree.size() == 15u );

                vector<uint64_t> expected_keys{0, 1, 1, 3, 2, 3, 3, 7, 4, 5, 5, 7, 6, 7, 7};
                for (uint i=0 ; i<15 ; i++) {
                    EXPECT( zrmt.tree[i].first == expected_keys[i] );
                    EXPECT( zrmt.tree[i].second == 0u );
                }
                cout << "  OK" << endl;
            }

            SECTION( "Update" ) {
                cout << "  Update" << endl;

                zrmt.update(3, 2);
                vector<uint64_t> expected_values1{0, 0, 0, 2, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0};
                for (uint idx=0 ; idx<zrmt.tree.size() ; idx++)
                    EXPECT(zrmt.tree[idx].second == expected_values1[idx]);

                zrmt.update(5, 3);
                vector<uint64_t> expected_values2{0, 0, 0, 2, 0, 2, 2, 3, 0, 3, 3, 3, 0, 0, 0};
                for (uint idx=0 ; idx<zrmt.tree.size() ; idx++)
                    EXPECT(zrmt.tree[idx].second == expected_values2[idx]);
                cout << "  OK" << endl;
            }

            SECTION( "0-Range Max" ) {
                cout << "  0-Range Max" << endl;
                
                zrmt.update(3, 2);
                zrmt.update(5, 3);

                EXPECT( zrmt.zero_range(1) == 0u );
                EXPECT( zrmt.zero_range(3) == 2u );
                EXPECT( zrmt.zero_range(4) == 2u );
                EXPECT( zrmt.zero_range(6) == 3u );

                EXPECT( zrmt.first_max_key(0) == 0u );
                EXPECT( zrmt.first_max_key(2) == 3u );
                EXPECT( zrmt.first_max_key(3) == 5u );

                cout << "  OK" << endl;
            }
        }

    }
};

extern lest::tests & specification();

MODULE( specification(), module )
