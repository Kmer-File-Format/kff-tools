// C++11 - use multiple source files.

#include "lest.hpp"
#include "RMT.hpp"
#include "compact.hpp"


using namespace std;

inline bool operator<(const PairInt& l, const PairInt& r)
{
    if (l.second == r.second)
        return l.first < r.first;
    else
        return l.second < r.second;
}

inline ostream& operator<<(ostream& os, const PairInt& p)
{
    os << p.first << ", " << p.second;
    return os;
}

const lest::test module[] = {

    CASE("Test a 0-Range Max Tree implementation") {
        cout << "Test 0-Range Max Tree" << endl;
        
        SETUP("Init") {
            vector<uint64_t> keys{0, 1, 2, 3, 4, 5, 6, 7};
            RangeMaxTree<uint64_t> rmt(keys);

            SECTION( "Test 2 Bytes bin" )
            {
                cout << "\tConstruction" << endl;
                EXPECT( rmt.tree.size() == 15u );

                vector<uint64_t> expected_keys{0, 1, 1, 3, 2, 3, 3, 7, 4, 5, 5, 7, 6, 7, 7};
                for (uint i=0 ; i<15 ; i++) {
                    EXPECT( rmt.tree[i].first == expected_keys[i] );
                    EXPECT( rmt.tree[i].second == 0u );
                }
                cout << "\t\tOK" << endl;
            }

            SECTION( "Update" ) {
                cout << "\tUpdate" << endl;

                rmt.update(3, 2);
                vector<uint64_t> expected_values1{0, 0, 0, 2, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0};
                for (uint idx=0 ; idx<rmt.tree.size() ; idx++)
                    EXPECT(rmt.tree[idx].second == expected_values1[idx]);

                rmt.update(5, 3);
                vector<uint64_t> expected_values2{0, 0, 0, 2, 0, 2, 2, 3, 0, 3, 3, 3, 0, 0, 0};
                for (uint idx=0 ; idx<rmt.tree.size() ; idx++)
                    EXPECT(rmt.tree[idx].second == expected_values2[idx]);
                cout << "\t\tOK" << endl;
            }

            SECTION( "0-Range Max" ) {
                cout << "\t0-Range Max" << endl;
                
                rmt.update(3, 2);
                rmt.update(5, 3);

                EXPECT( rmt.zero_range(1) == 0u );
                EXPECT( rmt.zero_range(3) == 2u );
                EXPECT( rmt.zero_range(4) == 2u );
                EXPECT( rmt.zero_range(6) == 3u );

                EXPECT( rmt.first_max_key(0) == 0u );
                EXPECT( rmt.first_max_key(2) == 3u );
                EXPECT( rmt.first_max_key(3) == 5u );

                cout << "\t\tOK" << endl;
            }

            SECTION( "Range Max" ) {
                cout << "\tRange Max" << endl;
                
                rmt.update(1, 7);
                rmt.update(2, 2);
                rmt.update(3, 9);
                rmt.update(5, 10);

                EXPECT( rmt.range(1, 2) == 7u );
                EXPECT( rmt.range(0, 2) == 7u );
                EXPECT( rmt.range(2, 2) == 2u );
                EXPECT( rmt.range(2, 4) == 9u );
                EXPECT( rmt.range(6, 7) == 0u );

                cout << "\t\tOK" << endl;
            }

            SECTION( "Bounded first max" ) {
                cout << "\tBounded first max" << endl;
                
                rmt.update(1, 7);
                rmt.update(2, 2);
                rmt.update(3, 9);
                rmt.update(5, 10);
                rmt.update(6, 7);

                EXPECT( rmt.bounded_first_max_key(9, 1) == 3u );
                EXPECT( rmt.bounded_first_max_key(10, 1) == 5u );
                EXPECT( rmt.bounded_first_max_key(7, 1) == 1u );
                EXPECT( rmt.bounded_first_max_key(7, 2) == 3u );

                cout << "\t\tOK" << endl;
            }
        }

        SETUP("With pair") {

            vector<PairInt> keys{PairInt(4, 0), PairInt(1, 1), PairInt(2, 1), PairInt(3, 1), PairInt(1, 3), PairInt(3, 3), PairInt(0, 5)};
            RangeMaxTree<PairInt> rmt(keys);
            uint64_t pos = 10;

            SECTION("find_second_coordinates") {
                cout << "\tfind_second_coordinates" << endl;
                EXPECT(rmt.find_second_coordinates(PairInt(4, 0), true) / 2u == 0u);
                EXPECT(rmt.find_second_coordinates(PairInt(4, 0), false) / 2u == 0u);

                EXPECT(rmt.find_second_coordinates(PairInt(1, 1), true) / 2u == 1u);
                EXPECT(rmt.find_second_coordinates(PairInt(1, 1), false) / 2u == 3u);

                EXPECT(rmt.find_second_coordinates(PairInt(1, 3), true) / 2u == 4u);
                EXPECT(rmt.find_second_coordinates(PairInt(1, 3), false) / 2u == 5u);

                EXPECT(rmt.find_second_coordinates(PairInt(2, 2), true) / 2u == 4u);
                EXPECT(rmt.find_second_coordinates(PairInt(2, 2), false) / 2u == 3u);

                EXPECT(rmt.find_second_coordinates(PairInt(2, 4), true) / 2u == 6u);
                EXPECT(rmt.find_second_coordinates(PairInt(2, 4), false) / 2u == 5u);

                cout << "\t\tOK" << endl;
            }

            SECTION("range_between") {
                cout << "\trange_between" << endl;

                rmt.update(PairInt(1, 1), 7u);
                rmt.update(PairInt(3, 1), 2u);
                rmt.update(PairInt(1, 3), 9u);
                rmt.update(PairInt(0, 5), 10u);

                EXPECT(rmt.range_between(PairInt(1, 1), PairInt(1, 1), pos) == 7u); EXPECT(pos == 1u);
                EXPECT(rmt.range_between(PairInt(1, 0), PairInt(1, 3), pos) == 9u); EXPECT(pos == 4u);
                EXPECT(rmt.range_between(PairInt(1, 3), PairInt(1, 3), pos) == 9u); EXPECT(pos == 4u);
                EXPECT(rmt.range_between(PairInt(1, 3), PairInt(0, 5), pos) == 10u); EXPECT(pos == 6u);
                EXPECT(rmt.range_between(PairInt(0, 5), PairInt(0, 5), pos) == 10u); EXPECT(pos == 6u);

                cout << "\t\tOK" << endl;
                }
            }
        }

    };

extern lest::tests & specification();

MODULE( specification(), module )
