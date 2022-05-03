

template <typename K, typename V>
class RMQ {
private:
		uint size;
    V nodes[];
    V leaves[];
		unordered_map<K, uint> positions;
public:
    RMQ(uint size);
    void update(K key, V val);
    V max(uint left, uint right);
};


template <typename V>
RMQ::RMQ(uint size) {
	this->size = size;
	
};

void RMQ::update(K, ) {

}
