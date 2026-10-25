/*
*****************************************************************************************************************
*****************************************************************************************************************
* De Bruijn Graph Genome Assembler
*
* A genome assembler that constructs a de Bruijn graph from sequencing
* reads, using read cleaning and graph simplification to reduce errors.
*
* Tip removal and bubble detection are used to identify and remove
* structures associated with sequencing errors.
*
* For paired-end reads, contig linking uses the estimated distance
* between read pairs to guide contig reconstruction.
*****************************************************************************************************************
*****************************************************************************************************************
*/

#include<vector>
#include<string>
#include<cstring>
#include<unordered_map>
#include<algorithm>
#include<unordered_set>
#include<utility>
#include<functional>
#include<queue>
#include<cstdint>
#include<iostream>

#include<fstream>


using namespace std;




struct STRING_REF{
    STRING_REF(const string*s,int st,int len):str(s),start(st),length(len){}
    STRING_REF() = default;
    const string* str;
    uint16_t start;
    uint16_t length;
    bool operator==(const STRING_REF&other)const{
        if(length != other.length){return false;}
        if(str == other.str && start == other.start){return true;}
        return(
            std::memcmp(str->data() + start, other.str->data() + other.start, length) == 0
        );
    }
    const char& operator[](size_t i)const{
        return(str->data()[start+i]);
    }
    void print_last_char()const{
        cout.write(str->data()+start+length-1,1);
    }
    friend ostream& operator<<(ostream&out,const STRING_REF& ref);
};


ostream& operator<<(ostream&out,const STRING_REF& ref){
    if(ref.str && ref.length>0){
        out.write(ref.str->data() + ref.start, ref.length);
    }
    return out;
}



class K_MER_BIT_MAP{
    public:
        K_MER_BIT_MAP(
            vector<string>&reads,
            vector<uint16_t>&dist,
            bool is_pair
        ){
            process_reads(reads,dist,is_pair);
        }

        K_MER_BIT_MAP() = default;

        typedef pair<uint64_t,uint64_t> K_MER_128;

        //Used to signal if a node or kmer is not found in the array
        static constexpr size_t NOT_FOUND = static_cast<size_t>(-1);

        struct pair_data{
            pair_data(uint32_t i,uint16_t d):target_index(i),dist(d){}
            uint32_t target_index;
            uint16_t dist;
        };
        unordered_map<int,vector<pair_data>> pair_map;



        size_t size()const{
            return arr.size();
        }

        string str_at(const size_t& i)const{
            if(i>=arr.size()){
                cout<<"size: "<<arr.size()<<" i: "<<i;
                throw std::out_of_range("Bit index out of range");
            }
            return decode_bit_to_str(i);
        }

        void print_str_at(const size_t& i,ostream&out)const{
            if(i>=arr.size()){
                cout<<"size: "<<arr.size()<<" i: "<<i;
                throw std::out_of_range("Bit index out of range");
            }
            print_bit_as_str(i,out);

        }

        K_MER_128 operator[](const size_t& i)const{
            return bit_at(i);
        }

        K_MER_128 bit_at(const size_t& i)const{
            if(i>=arr.size()){
                throw std::out_of_range("Bit index out of range");
            }
            return arr[i];
        }



        size_t find(const K_MER_128& k_mer)const{
            auto it = lower_bound(arr.begin(),arr.end(),k_mer);
            if(it != arr.end() && *it == k_mer){
                return (distance(arr.begin(),it));
            }
            else{
                return NOT_FOUND;
            }
        }

        size_t find(const STRING_REF&ref)const{
            K_MER_128 bit = encode_str_to_bit(ref);
            return find(bit);
        }

        size_t end()const{
            return NOT_FOUND;
        }



        void print_back_char_at(size_t i,ostream&out)const{
            if(i>=arr.size()){
                throw std::out_of_range("Bit index out of range");
            }
            //since I want the last char I
            //use the low side and I shift by 0
            //like in (00 10 11) I go to 11 
            //transform it to its base (T in this case) and print it
            out<<get_base(arr[i].second,0);
        }

        void print_mem_size(){
            cout<<"\nBIT_MAP_SIZE: "<<
            (
                sizeof(vector<K_MER_128>) + sizeof(K_MER_128) * arr.size()+
                sizeof(vector<uint8_t>) + sizeof(uint8_t) * count_arr.size()
            ) / (1024.0 * 1024) <<" MB\n";
        }

        uint16_t get_k_len(){
            return k_mer_len + 1;
        }


        struct K_MER_HASHER{
            size_t operator()(const K_MER_128& k_mer)const{
                size_t h1 = hash<uint64_t>{}(k_mer.first);
                size_t h2 = hash<uint64_t>{}(k_mer.second);
                return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
            }
        };


        void build_pair_link_map(
            const vector<string>&reads,
            const vector<uint16_t>&distances,
            const unordered_set<int>&junction_nodes,
            int offset = 1, bool check_set = true
        ){
            for(int i = 0; i+1<reads.size(); i+=2){
                if(reads[i].size() < k_mer_len || reads[i + 1].size() < k_mer_len){
                    continue;
                }
                int size = min(reads[i].size(),reads[i+1].size());
                for(int start_pos = 0; start_pos<= size-k_mer_len; start_pos+=offset){
                    K_MER_128 k1 = encode_str_to_bit(STRING_REF(&reads[i],start_pos,k_mer_len));
                    size_t k1_index = find(k1);
                    if(check_set && !junction_nodes.count(k1_index)){continue;}
                    K_MER_128 k2 = encode_str_to_bit(STRING_REF(&reads[i+1],start_pos,k_mer_len));
                    size_t k2_index = find(k2);
                    
                    if(k2_index != NOT_FOUND && k1_index != NOT_FOUND){
                        //cout<<"pair: "<<k1_index<<"->"<<k2_index<<" dist: "<<distances[i/2]<<"\n";
                        //the index at arr acts as the verex id, 
                        //i.e if the vertex is 5 then its at index 5 in the arr
                        pair_map[k1_index].push_back({static_cast<uint32_t>(k2_index),distances[i/2]});
                    }
                }
            }
        }




    private:

        uint16_t k_mer_len;

        //must be in this order because of how I encoded the bit
        string bases = "ACGT";

        vector<K_MER_128> arr;
        vector<uint8_t>count_arr;
        int16_t min_freq = 2;

        static const uint8_t big_k_mer = 50;
        static const uint8_t low_k_mer = 20;

        size_t gene_vol = 0;
        size_t unique_gene = 0;


        void reset(){
            arr.clear(); arr.shrink_to_fit();
            count_arr.clear(); count_arr.shrink_to_fit();
            gene_vol = 0;
            unique_gene = 0;
        }


        uint16_t count_at(size_t i){
            if(i>=arr.size()){
                throw std::out_of_range("Bit index out of range");
            }
            return count_arr[i];
        }

        uint16_t count_of_k_mer(const STRING_REF&ref){
            size_t i = find(ref);
            if(i != NOT_FOUND){
                return count_arr[i];
            }
            else{
                return 0;
            }
        }


        void process_reads(vector<string>&reads,vector<uint16_t>&dist,bool is_pair){
            k_mer_len = low_k_mer;
            clean_reads(reads,dist,is_pair);
            reset();
            k_mer_len = big_k_mer;
            create_bit_arr(reads,dist,is_pair,false);

            //create general link map for all the reads with an offset of 10
            //so not all kmers would be included to save memory
            unordered_set<int> empty_set;
            build_pair_link_map(reads,dist,empty_set,20,false);

            //to save memory I can clear the dist
            //dist.clear(); dist.shrink_to_fit();
        }



        size_t clean_reads(
            vector<string>&reads,
            vector<uint16_t>&dist,
            bool is_pair
        ){
            size_t reads_cleaned = 0;
            bool cleaning_reads = true;
            
            create_bit_arr(reads,dist,is_pair,cleaning_reads);
            for(string&e:reads){
                K_MER_128 current_k_mer = encode_str_to_bit(STRING_REF(&e,0,k_mer_len));
                bool first_is_error = (find(current_k_mer) == NOT_FOUND);
                for(int i = k_mer_len; i<=e.size()-k_mer_len; i++){
                    current_k_mer = append_char_to_k_mer(current_k_mer,e[i]);
                    bool is_last_char = true;
                    if(!first_is_error){
                        if(find(current_k_mer) == NOT_FOUND){
                            char base = find_best_char_for_erroneous_kmer(current_k_mer,is_last_char);
                            if(base != '\0'){
                                e[i] = base;
                                reads_cleaned++;
                                break;
                            }
                            
                        }
                    }
                    else{
                        if(find(current_k_mer) != NOT_FOUND){
                            is_last_char = false;
                            //add a push to front of kmer function for better speed
                            current_k_mer = push_char_to_front_k_mer(current_k_mer,e[i-k_mer_len]);
                            char base = find_best_char_for_erroneous_kmer(current_k_mer,is_last_char);
                            if(base != '\0'){
                                e[i-k_mer_len] = base;
                                reads_cleaned++;
                            }
                            break;
                        }
                    }
                }
            }
            //cout<<"READS CLEANED = "<<reads_cleaned<<endl;
            return reads_cleaned;
        }

        char find_best_char_for_erroneous_kmer(K_MER_128& k_mer,bool is_last_char){
            uint16_t best_count = 0;
            uint64_t high_len = k_mer_len/2;
            uint64_t low_len = k_mer_len - high_len;
            char best_char = '\0';
            for(uint64_t bit_base = 0; bit_base<4; bit_base++){
                if(is_last_char){
                    k_mer.second = (k_mer.second & ~3ULL) | bit_base;
                }
                else{
                    uint64_t shift = (high_len - 1)*2;
                    k_mer.first = (k_mer.first & ~(3ULL<<shift)) | (bit_base<<shift);
                }

                size_t index = find(k_mer);
                if(index != NOT_FOUND){
                    uint16_t count = count_at(index);
                    if(best_count < count){
                        best_count = count;
                        best_char = bases[bit_base];
                    }
                }
            }
            return best_char;
        }

        void create_bit_arr(
            vector<string>&reads,
            vector<uint16_t>&dist,
            bool is_pair,bool cleaning_reads
        ){
            uint64_t high_len = k_mer_len/2;
            uint64_t prefix_shift = (2 * high_len) - 4;
            for(uint64_t prefix = 0; prefix<16; prefix++){
                vector<K_MER_128> chunk;
                chunk.reserve(7000000);
                for(const string& e:reads){
                    if(e.size()<k_mer_len){
                        continue;
                    }
                    K_MER_128 current_k_mer = encode_str_to_bit(STRING_REF(&e,0,k_mer_len));
                    if((current_k_mer.first>>prefix_shift) == prefix){
                        chunk.push_back(current_k_mer);
                    }
                    for(size_t i = k_mer_len; i<e.size(); i++){
                        current_k_mer = append_char_to_k_mer(current_k_mer,e[i]);
                        if((current_k_mer.first>>prefix_shift) == prefix){
                            chunk.push_back(current_k_mer);
                        }
                    }
                }
                
                sort_k_mer(chunk);
                clean_bit_arr(chunk,cleaning_reads);

                if(prefix == 0 && !cleaning_reads){
                    bool restart_needed = tune_settings();
                    if(restart_needed){
                        reset();
                        create_bit_arr(reads,dist,is_pair,false);
                        return;
                    }
                }

                arr.insert(arr.end(),chunk.begin(),chunk.end());
            }
        }


        

        K_MER_128 encode_str_to_bit(const STRING_REF&ref)const{
            uint64_t low = 0;
            uint64_t high = 0;
            uint16_t half_len = k_mer_len/2;
            for(uint16_t i = 0; i<k_mer_len; i++){
                uint64_t base = get_bit_id(ref[i]);
                
                //this shifts the high/low by 2 bits to the left
                //adding 00 to the right which we populate with base
                //so if we start with G then add T it'll be like this:
                //10 -> 1000 -> 1011
                if(i < half_len){
                    high = (high << 2) | base;
                }
                else{
                    low = (low << 2) | base;
                }
            }
            return{high,low};
        }

        string decode_bit_to_str(const size_t& k_mer_indx)const{
            uint16_t half_len = k_mer_len / 2;
            uint16_t second_half_len = k_mer_len - half_len;
            string result(k_mer_len,'\0');

            for(uint16_t i = 0; i<half_len; i++){
                //the shift needed is how many bits away
                //I am from the far right. this is calculated by 
                //how many bits I shifted it to get to that pos beforehand
                //so for 00 10 01 11 lets say I want 10
                //the shift is 4 -> 2 *( half_len (4) - 1 - i (1) )
                uint16_t shift = 2 * (half_len-1-i);
                result[i] = get_base(arr[k_mer_indx].first, shift);
            }
            for(uint16_t i = 0; i<second_half_len; i++){
                //same thing as before just for the low
                //part of the pair instead of the high
                uint16_t shift = 2 * (second_half_len-1-i);
                result[half_len+i] = get_base(arr[k_mer_indx].second, shift);
            }
            return result;
        }

        //same logic as the decode to str but Here I just print it
        void print_bit_as_str(const size_t& k_mer_indx,ostream&out)const{
            uint16_t half_len = k_mer_len / 2;
            uint16_t second_half_len = k_mer_len - half_len;

            for(uint16_t i = 0; i<half_len; i++){
                uint16_t shift = 2 * (half_len-1-i);
                out<<get_base(arr[k_mer_indx].first, shift);
            }
            for(uint16_t i = 0; i<second_half_len; i++){
                uint16_t shift = 2 * (second_half_len-1-i);
                out<<get_base(arr[k_mer_indx].second, shift);
            }
        }


        char get_base(const uint64_t& half_k_mer,const uint16_t& shift = 0)const{
            //shift that bit by shift amount to the right
            //making it at the far right where I get rid of everything
            //above the first 2 bits i.e every thing bigger than 3
            //giving me the index of the base from what I encoded it
            //i.e 10 = 2 = G which matches the bases' order
            uint16_t base_index = (half_k_mer >> shift) & 3;
            return bases[base_index];
        }

        K_MER_128 push_char_to_front_k_mer(K_MER_128 k_mer,const char& c){
            uint16_t high_len = k_mer_len / 2;
            uint16_t low_len = k_mer_len - high_len;

            uint64_t base = get_bit_id(c);

            //because I will push to the front I need to get
            //the last char in the high to move it to the front of low
            uint64_t shifted_from_high = k_mer.first & 3ULL;
            //pushes the char to the front of high while removing the last 2 bits
            k_mer.first = (base<<((high_len-1 )*2)) | (k_mer.first >> 2);
            //adding the last two bits of high to front of low
            k_mer.second = (shifted_from_high<<((low_len-1)*2)) | (k_mer.second >> 2);
            return k_mer;
        }

        K_MER_128 append_char_to_k_mer(K_MER_128 k_mer,const char& c){
            uint16_t high_len = k_mer_len / 2;
            uint16_t low_len = k_mer_len - high_len;

            //creates 00111...111 so when I do & mask the right most will be removed
            uint64_t high_mask = (1ULL << (2*high_len)) -1;
            uint64_t low_mask = (1ULL << (2*low_len)) -1;

            uint64_t base = get_bit_id(c);

            //gets the first two bits that now because
            //of appending to low need to be moved to back of high
            uint64_t shifted_from_low = (k_mer.second >> (2*(low_len-1))) & 3ULL;
            //appends the new base to the left of low and removes the 2 extra at the right
            k_mer.second = ((k_mer.second << 2) | base) & low_mask;
            //appends the 2 extra from low to the left of high and removes the two exta at the right
            k_mer.first = ((k_mer.first << 2) | shifted_from_low) & high_mask;
            return k_mer;
        }

        uint64_t get_bit_id(const char& c)const{
            uint64_t id = 0;
            switch(c){
                case 'A': id = 0; break;
                case 'C': id = 1; break;
                case 'G': id = 2; break;
                case 'T': id = 3; break;
                default: id = 0;break;
            }
            return id;
        }

        void sort_k_mer(vector<K_MER_128>&a){
            sort(a.begin(),a.end());
        }

        double coverage_of_k(){
            return ( gene_vol / static_cast<double>(unique_gene) );
        }

        bool tune_settings(){

            double CK = coverage_of_k();

            double target_ck = k_mer_len > 30 ? 12.0 : 6.0;
            if( CK >= target_ck){
                if(k_mer_len % 2 == 1){
                    k_mer_len -= 1;
                }
                else{
                    //already good coverage no need to change anything
                    return false;
                }
                
            }
            else if(k_mer_len > low_k_mer){
                min_freq = 2;
                int next_k_len = k_mer_len - 5;
                int min_k = low_k_mer;
                k_mer_len = max(min_k,next_k_len);
            }
            else if(min_freq != 1){
                min_freq = 1;
                k_mer_len = low_k_mer;
            }
            else{
                //ck is low but the best we can do so return false;
                return false;
            }
            //if(k_mer_len%2 == 1){k_mer_len -= 1;}
            //cout<<"CHANGE NEEDED - chosen settings are: K len = "<<k_mer_len<<" min freq = "<<min_freq<<"\n";
            return true;
        }

        

        void clean_bit_arr(vector<K_MER_128>&a,bool is_cleaning_reads = false){
            
            size_t write_index = 0;
            size_t size = a.size();
            for(size_t i = 0; i<size; ){
                size_t j = i+1;
                while(j<size && a[i] == a[j]){
                    j++;
                }
                uint16_t freq = j-i;

                if(freq >= 2){
                    gene_vol += freq;
                    unique_gene++;
                }

                if(freq >= min_freq){
                    a[write_index] = a[i];
                    write_index++;
                    if(is_cleaning_reads){
                        count_arr.push_back(freq);
                    }
                    
                }
                i = j;
            }
            a.resize(write_index);
            a.shrink_to_fit();
        }
        

};



struct edge{
    edge(int t,int w):to(t),weight(w),visited(false){}
    edge() = default;
    int to;
    unsigned short weight;
    bool visited = false;
    
    void operator++(){
        weight++;
    }
};

struct DE_BRUIJN_ROW{
    DE_BRUIJN_ROW(vector<edge>r):row(r){}
    DE_BRUIJN_ROW() = default;
    vector<edge>row;
    edge& operator[](int i){
        if(i>=row.size()){throw std::out_of_range("Column index out of range");}
        return row[i];
    }
    size_t size()const{
        return row.size();
    }
    void push_back_edge(const edge& e){
        row.push_back(e);
    }
    void erase(const int& i){
        row.erase(row.begin()+i);
    }
    bool empty()const{
        return row.empty();
    }
};


template<typename K, typename V, typename H>
void print_unordered_map_mem_size(const unordered_map<K,V,H>&map){
    int node_overhead = 32;
    cout<<"\nmap size: "<<
    (
        sizeof(unordered_map<K,V,H>) + (map.bucket_count()*sizeof(void*)) +
        (map.size() * (sizeof(K) + sizeof(V) + node_overhead)) + sizeof(H)
    )/(1024.0*1024)
    <<" MB\n";
}



class DE_BRUIJN_GRAPH{
    public:

        DE_BRUIJN_GRAPH(vector<string>& entries,vector<uint16_t>&dist,bool is_pair){
            create_k_mer_graph(entries,dist,is_pair);
            
        }
        DE_BRUIJN_GRAPH() = default;
        
        vector<DE_BRUIJN_ROW>graph;
        
        vector<int>in_deg,out_deg;
        
        K_MER_BIT_MAP id_to_str;
        uint16_t true_k_len;
        
        DE_BRUIJN_ROW& operator[](int v){
            if(v>=graph.size()){throw std::out_of_range("Row index out of range");}
            return graph[v];
        }
        
        
        void update_edge_degree(){
            in_deg.assign(graph.size(),0);
            out_deg.assign(graph.size(),0);
            size_t size = graph.size();
            for(size_t v = 0; v<size; v++){
                for(size_t j = 0; j<graph[v].size(); j++){
                    size_t u = graph[v][j].to;
                    in_deg[u]++;
                    out_deg[v]++;
                }
            }
        }
        
        
        int find_edge(const int& v,const int& u){
            for(int i = 0; i<graph[v].size(); i++){
                if(graph[v][i].to == u){
                    return i;
                }
            }
            return -1;
        }
        
        void remove_edge(const int& v, const int& u){
            int i = find_edge(v,u);
            if(i==-1){return;}
            graph[v].erase(i);
            in_deg[u]--;
            out_deg[v]--;
        }

        void remove_edge_by_index(int v, int i){
            int u = graph[v][i].to;
            graph[v].erase(i);
            in_deg[u]--;
            out_deg[v]--;
        }
        
        size_t size()const{
            return graph.size();
        }

        size_t get_total_edges(){
            size_t total = 0;
            for(int i = 0; i<graph.size(); i++){
                total += graph[i].size();
            }
            return total;
        }

        void clear(){
            graph.clear();
            graph.shrink_to_fit();
        }


        //removes Tips and bubbles over and over again until none where removed
        void remove_tips_and_bubbles(int max_depth){

            TIP_REMOVER tip_remover(*this,max_depth);

            BUBBLE_REMOVER bubble_remover(*this,max_depth);

            unsigned int total = 0;

            while(true){
                unsigned int tips = 0;
                unsigned int bubbles = 0;
                tips += tip_remover.remove_tips();
                bubbles += bubble_remover.remove_bubbles();
                //cout<<"Total tips this round: "<<tips<<" Total bubbles this round: "<<bubbles<<"\n";
                total += tips+bubbles;
                if(tips+bubbles == 0){
                    break;
                }
            }
            //cout<<"Total in all rounds: "<<total<<"\n";
            
        }


        void remove_relative_low_weight_edges(){
            //remove an edge if its weight extreamly low than the rest of the edges
            for(int v = 0; v<graph.size(); v++){
                int max_weight = 0;
                for(int i = 0; i<graph[v].size(); i++){
                    max_weight =(max_weight < graph[v][i].weight ? graph[v][i].weight : max_weight);
                }
                for(int i = graph[v].size()-1; i>=0; i--){
                    if(graph[v][i].weight < 0.2 * max_weight){
                        remove_edge_by_index(v,i);
                    }
                }
            }  
        }


        void build_pair_link_map(const vector<string>&reads,const vector<uint16_t>&dist){
            unordered_set<int> junction_nodes;
            for(int v = 0; v<graph.size(); v++){
                if(graph[v].size()>=2){
                    junction_nodes.insert(v);
                }
            }
            id_to_str.build_pair_link_map(reads,dist,junction_nodes);
        }
    
        
        void print_graph(){
            for(int v = 0; v<graph.size(); v++){
                for(int i = 0; i<graph[v].size(); i++){
                    cout<<v<<"->"<<graph[v][i].to<<" = "
                    <<id_to_str.str_at(v)<<"->"<<id_to_str.str_at(graph[v][i].to)
                    <<" Weight: "<<graph[v][i].weight<<"\n";
                }
            }
        }
        

        void print_graph_mem_size(){
            size_t total = 0;
            for(int i = 0; i<graph.size(); i++){
                total += get_row_mem_size(graph[i]);
            }
            cout<<"\ngraph size: "<<total/(1024.0*1024)<<" MB\n";
        }

    
    private:
    
        void create_k_mer_graph(
            vector<string>& entries,
            vector<uint16_t>&dist,
            bool is_pair
        ){
            int count = entries.size();

            //create an the kmer map to be used to build the graph
            id_to_str = K_MER_BIT_MAP(entries,dist,is_pair);

            true_k_len = id_to_str.get_k_len();

            size_t size = id_to_str.size();
            graph.resize(size);


            for(const string& e:entries){
                if(e.size()<true_k_len){continue;}
                //sliding a window through the entry to get all k-mers
                for(int i = 0; i<=e.size()-true_k_len; i++){

                    //vert size is k-1 so the edge is k long
                    STRING_REF pre(&e,i,true_k_len-1);
                    STRING_REF suff(&e,i+1,true_k_len-1);
                    
                    size_t u = id_to_str.find(pre);
                    size_t v = id_to_str.find(suff);
                    if(u == id_to_str.end()){continue;}
                    if(v == id_to_str.end()){continue;}


                    //prevent the same edge from being created twice
                    int v_index = find_edge(u,v);
                    if(v_index == -1){
                        add_edge(u,v);
                        v_index = graph[u].size() - 1;
                    }
                    
                    //++ increments the weight of the edge
                    ++graph[u][v_index];
                }
            }
        }

        
        void add_row(const DE_BRUIJN_ROW& row){
            graph.push_back(row);
        }
        void add_row(){
            DE_BRUIJN_ROW row;
            graph.push_back(row);
        }
        
        
        //returns edge from a specific v to a specific u
        edge* get_edge(const int& v,const int& u){
            int i = find_edge(v,u);
            if(i != -1){
                return &graph[v][i];
            }
            throw std::out_of_range("Edge not found!");
            
            return nullptr;
        }
        
        void add_edge(const int&v,const int u){
            edge e(u,0);
            graph[v].push_back_edge(e);
        }
        

        size_t get_row_mem_size(DE_BRUIJN_ROW& row){
            const vector<edge>& vec = row.row;
            return(sizeof(vector<edge>) + sizeof(edge)*vec.size());
        }


        
        //DE BRUIJN GRAPH's private classes


//##############################################################################################################
//######################################### BUBBLE REMOVER START ###############################################
//##############################################################################################################   

        
        class BUBBLE_REMOVER{
            public:
            BUBBLE_REMOVER(DE_BRUIJN_GRAPH&g,int depth):
            graph(g),in_deg(g.in_deg),out_deg(g.out_deg),max_depth(depth){}
            
            ~BUBBLE_REMOVER() = default;
            
            
            int remove_bubbles(){
                graph.update_edge_degree();
                unordered_set<int>in_cannidates,out_cannidates;
                record_valid_bubble_vertices_cannidates(in_cannidates,out_cannidates);
                return pop_valid_bubbles(in_cannidates,out_cannidates);
            }
            
            
            private:
            DE_BRUIJN_GRAPH& graph;
            vector<int>&in_deg;
            vector<int>&out_deg;
            int max_depth;
            
            
            void record_valid_bubble_vertices_cannidates(
                unordered_set<int>&in_cannidates,unordered_set<int>&out_cannidates
            ){
                int size = in_deg.size();
                for(int v = 0; v<size; v++){
                    if(in_deg[v]>=2){
                        in_cannidates.insert(v);
                    }
                    if(out_deg[v]>=2){
                        out_cannidates.insert(v);
                    }
                }
            }
            
            bool disjoint_paths(const vector<int>&a,const vector<int>&b){
                vector<bool>visited(graph.size(),false);
                int a_size = a.size();
                int b_size = b.size();
                for(const int& u:a){
                    visited[u] = true;
                }
                //remark the shared v and w as false
                visited[a[0]] = visited[a.back()] = false;
                for(const int u:b){
                    if(visited[u]){return false;}
                }
                return true;
            }
            
            void find_paths(
                const int& curr,
                vector<int>&path,
                vector<bool>&visited,
                unordered_set<int>&in_cannidates,
                unordered_map<int,vector<vector<int>>>& to_w_paths
            ){
                //depth is in edges count not vertex count
                int depth = path.size()-1;
                auto it = in_cannidates.find(curr);
                if(depth>0 &&  it != in_cannidates.end()){
                    to_w_paths[curr].push_back(path);
                    //we dont return here because we might hit another w later
                    //before reaching the depth limit
                }
                
                if(depth == max_depth){return;}
                
                for(int i = 0; i<graph[curr].size(); i++){
                    const int& u = graph[curr][i].to;
                    if(visited[u]){continue;}
                    visited[u] = true;
                    path.push_back(u);
                    
                    find_paths(u,path,visited,in_cannidates,to_w_paths);
                    
                    //we track back here so we could use this
                    //vertex for another path to some other w
                    path.pop_back();
                    visited[u] = false;
                }
            }
            
            double get_path_weight(const vector<int>&path){
                int size = path.size();
                int edge_count = size-1;
                double total_weight = 0.0;
                for(int i = 0; i<size-1; i++){
                    const int& v = path[i];
                    const int& u = path[i+1];
                    int j = graph.find_edge(v,u);
                    if(j != -1){
                        total_weight += graph[v][j].weight;
                    }
                }
                return total_weight / edge_count;
            }
            void remove_path(const vector<int>&path){
                int size = path.size();
                for(int i = 0; i< size-1; i++){
                    const int& u = path[i];
                    const int& v = path[i+1];
                    graph.remove_edge(u,v);
                }
            }
            
            int pop_valid_bubbles(
                unordered_set<int>&in_cannidates,unordered_set<int>&out_cannidates
            ){
                int bubble_count = 0;
                int vert_count = graph.size();
                vector<vector<int>> paths_to_remove;
                
                for(const int& v:out_cannidates){
                    //keeps track af all paths from v to some w (an in_candidate)
                    unordered_map<int,vector<vector<int>>> to_w_paths;
                    
                    vector<bool>visited(vert_count,false);
                    visited[v] = true;
                    vector<int>path = {v};

                    
                    
                    //records all paths from v to w in the to_w_paths
                    find_paths(v,path,visited,in_cannidates,to_w_paths);
                    
                    //for each pair of paths from v to some w we check
                    //if all the inward nodes are different from each other
                    //(if the two paths are disjoint) if they are then a bubble exist
                    //so we count it
                    for(const auto& w_paths:to_w_paths){
                        //w_paths.second is the paths from v to w
                        int path_count = w_paths.second.size();
                        
                        for(int i = 0; i<path_count; i++){
                            for(int j = i+1; j<path_count; j++){
                                if(disjoint_paths(w_paths.second[i],w_paths.second[j])){
                                    double weight_i = get_path_weight(w_paths.second[i]);
                                    double weight_j = get_path_weight(w_paths.second[j]);
                                    const vector<int>&path_to_remove = (
                                    weight_i > weight_j ? w_paths.second[j] : w_paths.second[i]);
                                    
                                    /*
                                    cout<<"bubble path removed:\n";
                                    for(const int&v:path_to_remove){
                                    cout<<v<<" ";
                                    }
                                    cout<<"\n";
                                    //*/
                                    paths_to_remove.push_back(path_to_remove);
                                    bubble_count++;
                                }
                            }
                        }
                    }
                }
                for(const vector<int>& path:paths_to_remove){
                    remove_path(path);
                }

               return bubble_count;     
            }
                
                
        };          
            
            

//##############################################################################################################
//########################################## TIP REMOVER START #################################################
//##############################################################################################################   
            
        class TIP_REMOVER{
                public:
                TIP_REMOVER(DE_BRUIJN_GRAPH&g,int max_tip):
                graph(g),in_deg(g.in_deg),out_deg(g.out_deg),max_tip_size(max_tip){}
                
                int remove_tips(){
                    graph.update_edge_degree();
                    int forward_removed = remove_forward_tips();
                    int backward_removed = remove_backward_tips();
                    return backward_removed + forward_removed;
                }
                
                private:
                DE_BRUIJN_GRAPH& graph;
                vector<int>&in_deg;
                vector<int>&out_deg;
                int max_tip_size;
                
                
                //backward tips are tips that have 0 out degree
                //therefore we need to trace them from a node with
                //at least 2 out - this would mean it might lead to
                //a dead-end tip
                void find_backward_tip_origin(queue<int>&tips){
                    int vert_count = graph.size();
                    for(int v = 0; v<vert_count; v++){
                        if(out_deg[v]>=2){
                            tips.push(v);
                        }
                    }
                }
                
                //forward tips are tips that start with a node
                //that has 0 in degree. Not all of these will be
                //tips but they are the place to start tracing
                void find_forward_tip_origin(queue<int>&tips){
                    int vert_count = graph.size();
                    for(int v = 0; v<vert_count; v++){
                        if(in_deg[v]==0){
                            tips.push(v);
                        }
                    }
                }
                
                int remove_forward_tips(){
                    queue<int>tips;
                    find_forward_tip_origin(tips);
                    int tips_removed_this_round = 0;
                    while(!tips.empty()){
                        int tip = tips.front();
                        tips.pop();
                        tips_removed_this_round += forward_tip_removal(tip,tips);
                    }
                    return tips_removed_this_round;
                }
                
                int remove_backward_tips(){
                    int tips_removed_this_round = 0;
                    queue<int>tips;
                    find_backward_tip_origin(tips);
                    while(!tips.empty()){
                        int tip = tips.front();
                        tips.pop();
                        tips_removed_this_round += backward_tip_removal(tip);
                    }
                    return tips_removed_this_round;
                }
                
                
                
                int forward_tip_removal(const int& tip_start,queue<int>&tips){
                    //since the graph changes I need to make sure
                    //that what I had in the queue is still valid
                    //meaning the forward tip indeed has an in degree
                    //that is 0 and out that is 1
                    if(in_deg[tip_start] != 0 || out_deg[tip_start] != 1){return 0;}
                    
                    int curr = graph[tip_start][0].to;
                    vector<int>path = {tip_start,curr};
                    
                    build_tip_path(path,curr);
                    
                    const int&tip_end = path.back();
                    
                    int edge_count = path.size()-1;
                    int tips_removed = 0;
                    if(edge_count>0 && (in_deg[tip_end] > 1 || out_deg[tip_end] > 1 || out_deg[tip_end] == 0)){
                        
                        remove_tip(path);
                        tips_removed+=path.size()-1;
                        //add new tip to queue if the removal of the current one made another
                        if(in_deg[tip_end] == 0 && out_deg[tip_end] > 0){
                            tips.push(tip_end);
                        }
                        
                    }
                    return tips_removed;
                }
                
                int backward_tip_removal(const int& tip_origin){
                    int tips_removed = 0;
                    //since the graph changes I need to make sure
                    //that what I had in the queue is still valid
                    //meaning the origin still has out >= 2
                    if(out_deg[tip_origin] < 2){return tips_removed;}
                    
                    //need to look at all the edges that could be tips
                    for(int i = graph[tip_origin].size()-1; i>=0; i--){
                        
                        int curr = graph[tip_origin][i].to;
                        if(in_deg[curr] != 1){continue;}
                        vector<int> path = {tip_origin,curr};
                        
                        build_tip_path(path,curr);
                        
                        bool removed_tip = false;
                        
                        const int&tip_end = path.back();
                        if(out_deg[tip_end] == 0 && in_deg[tip_end] == 1){
                            
                            remove_tip(path);
                            
                            tips_removed+=path.size()-1;
                        }
                        
                    }
                    return tips_removed;
                }
                
                
                void build_tip_path(vector<int>&path,int curr){
                    /*
                    cout<<"curr path: "<<path[0]<<" "<<curr<<" ";
                    //*/
                    while(out_deg[curr] == 1 && in_deg[curr] == 1 && path.size() <= max_tip_size){
                        curr = graph[curr][0].to;
                        /*
                        cout<<curr<<" ";
                        //*/
                        path.push_back(curr);
                    }
                    /*
                    cout<<"\n";
                    //*/
                }
                
                void remove_tip(const vector<int>&path){
                    int edge_count = path.size()-1;
                    for(int i = 0; i<edge_count; i++){
                        const int& u = path[i];
                        const int& v = path[i+1];
                        graph.remove_edge(u,v);
                    }
                }
        };
            

//##############################################################################################################
//######################################## DE BRUIJN GRAPH CLASS END ###########################################
//##############################################################################################################  
    
};





    
class GENOME_ASSEMBLER{
    public:
        GENOME_ASSEMBLER(vector<string>&r,vector<uint16_t>&dist,bool is_pair)
        :reads(r),distances(dist),is_read_pair(is_pair){
        
            graph = DE_BRUIJN_GRAPH(reads,distances,is_read_pair);

            graph.remove_tips_and_bubbles(graph.true_k_len);
            graph.remove_relative_low_weight_edges();

            if(is_pair){
                graph.build_pair_link_map(reads,dist);
            }
            //can remove the reads as I don't need them anymore
            reads.clear(); reads.shrink_to_fit();
            dist.clear(); dist.shrink_to_fit();
        }
        
        void assemble_genome(){
            graph.update_edge_degree();
            const string contig_file_path = "genome_assembly/contig_output.txt";
            ofstream file(contig_file_path);
            build_contigs();

            if(is_read_pair){
                CONTIG_LINKER linker(contigs,graph);
                linker.join_linked_contigs();
            }

            print_contigs(file);
            cout<<"\nGenome assembled!\ncontigs outputed at "<<contig_file_path<<endl;
            
        }

        
    private:
        
        DE_BRUIJN_GRAPH graph;

        bool is_read_pair;
        vector<uint16_t>&distances;
        vector<string>&reads;


        struct CONTIG{
            CONTIG(vector<int>a):nodes(a){}
            vector<int>nodes;
            void insert(const CONTIG&other){
                nodes.insert(nodes.end(),other.nodes.begin(),other.nodes.end());
            }
            int& operator[](int i){
                if(i > nodes.size()){
                    throw std::out_of_range("contig base index out of range");
                }
                return nodes[i];
            }
        };
        vector<CONTIG>contigs;


        //start node that has out > 0 and is not a 1 to 1 node
        void find_start_edge(int& u,int& j){
            for(int v = 0; v<graph.size(); v++){
                if(graph.out_deg[v] > 0 && !(graph.out_deg[v] == 1 && graph.in_deg[v] == 1)){
                    for(int i = 0; i<graph[v].size(); i++){
                        if(!graph[v][i].visited){
                            u = v; j = i;
                            return;
                        }
                    }
                }
            }
        }

        int find_active_edge(int v){
            for(int i = 0; i<graph[v].size(); i++){
                if(!graph[v][i].visited){
                    return i;
                }
            }
            return -1;
        }


        
        void build_contigs(){

            //cout<<"building contigs...\n";

            //linear contigs starting from a start node that 
            //has out > 0 and !(in == 1 && out == 1)
            const K_MER_BIT_MAP&id_to_str = graph.id_to_str;

            
            const int min_contig_len = 1;

            while(true){
                vector<int> contig_path;
                int v = -1,u_index = -1;
                find_start_edge(v,u_index);
                if(v == -1){break;}
                const int&u = graph[v][u_index].to;
                graph[v][u_index].visited = true;
                contig_path.push_back(v);
                
                int curr = u;

                while(!graph[curr].empty()){
                    int next;
                    int edge_indx = -1;

                    contig_path.push_back(curr);

                    if(graph.in_deg[curr] == 1 && graph.out_deg[curr] == 1){
                        if(!graph[curr][0].visited){
                            edge_indx = 0;
                        }
                        else{break;}
                    }

                    if(edge_indx == -1){break;}

                    next = graph[curr][edge_indx].to;
                    graph[curr][edge_indx].visited = true;
                    /*
                    cout<<"->"<<next;
                    //*/
                    curr = next;
                }
                if(contig_path.size()>min_contig_len){
                    contigs.push_back(CONTIG(contig_path));
                }
                
                
            }
            
            //independent cyclic contigs
            for(int v = 0; v<graph.size(); v++){
                for(int i = 0; i<graph[v].size(); i++){
                    vector<int>contig_path;
                    if( graph[v][i].visited ){continue;}
                    
                    graph[v][i].visited = true;
                    contig_path.push_back(v);
                    
                    const int& target = graph[v][i].to;
                    
                    int curr = target;
                    
                    while(curr != v && !graph[curr].empty()){
                        contig_path.push_back(curr);
                        int edge_indx = find_active_edge(curr);
                        if(edge_indx == -1){break;}
                        const int& next = graph[curr][edge_indx].to;
                        graph[curr][edge_indx].visited = true;

                        curr = next;
                    }
                    if(contig_path.size() > min_contig_len){
                        contigs.push_back(CONTIG(contig_path));
                    }
                }
            }
        }


        void print_contigs(ostream&out){
            const K_MER_BIT_MAP&id_to_str = graph.id_to_str;
            size_t curr_contig = 0;
            for(int i = 0; i<contigs.size(); i++){
                //if(is_read_pair && joined_contigs[i]){continue;}
                curr_contig++;
                if(contigs[i].nodes.empty()){continue;}
                out<<">CONTIG"<<curr_contig<<"\n";
                id_to_str.print_str_at(contigs[i][0],out);
                for(int j = 1; j<contigs[i].nodes.size(); j++){
                    id_to_str.print_back_char_at(contigs[i][j],out);
                }
                out<<"\n";
            }
        }


        class CONTIG_LINKER{
            public:
                CONTIG_LINKER(vector<GENOME_ASSEMBLER::CONTIG>&cons,DE_BRUIJN_GRAPH& g):
                contigs(cons),graph(g)
                {
                    joined_contigs.assign(contigs.size(),false);
                    contig_link_matrix.resize(contigs.size(),vector<CONTIG_EDGE>(contigs.size()));
                }

                void join_linked_contigs(){
                    assign_node_to_contig_index();
                    build_contig_link_matrix();
                    //node to contig is not needed anymore so I can clear it
                    node_to_contig.clear(); node_to_contig.shrink_to_fit();
                    for(int i = 0; i<contig_link_matrix.size(); i++){
                        if(joined_contigs[i]){continue;}
                        int c = find_best_contig_match(i);
                        if(c != -1){
                            join_contigs(i,c);
                        }
                        
                    }
                    clear_joined_contigs();
                }
            
            private:
                vector<bool>joined_contigs;
                vector<GENOME_ASSEMBLER::CONTIG>&contigs;
                DE_BRUIJN_GRAPH&graph;
                struct NODE_CONTIG_DATA{
                    NODE_CONTIG_DATA(int i, int d):contig_index(i),depth(d){}
                    int contig_index;
                    int depth;
                };
                struct CONTIG_EDGE{
                    CONTIG_EDGE() = default;
                    int weight = 0;
                    int total_dist = 0;
                    int dist(){
                        return total_dist / weight;
                    }
                };
                vector<vector<NODE_CONTIG_DATA>>node_to_contig;
                vector<vector<CONTIG_EDGE>>contig_link_matrix;

                vector<int>contig_link;

                const int MIN_WEIGHT = 1;
                const int MAX_VARIANCE = 50;

                void assign_node_to_contig_index(){
                    node_to_contig.assign(graph.size(),{});
                    for(int i = 0; i<contigs.size(); i++){
                        for(int j = 0; j<contigs[i].nodes.size(); j++){
                            int v = contigs[i].nodes[j];
                            node_to_contig[v].push_back({i,j});
                        }
                    }
                }


                void build_contig_link_matrix(){
                    //cout<<"building links\n";
                    
                    for(auto it:graph.id_to_str.pair_map){
                        for(int i = 0; i<node_to_contig[it.first].size(); i++){
                            int contig_a_indx = node_to_contig[it.first][i].contig_index;
                            if(contig_a_indx>=0){

                                int contig_a_size = contigs[contig_a_indx].nodes.size();
                                int a_node_depth = node_to_contig[it.first][i].depth;

                                int a_dist_from_end =  contig_a_size - a_node_depth; 

                                for(const auto& data:it.second){
                                    for(int j = 0; j<node_to_contig[data.target_index].size(); j++){
                                        int contig_b_indx = node_to_contig[data.target_index][j].contig_index;
                                        if(contig_b_indx >= 0){
                                            int b_node_depth = node_to_contig[data.target_index][j].depth;
                                            int estimated_dist = data.dist;
                                            int bases_dist = (a_dist_from_end + b_node_depth) + graph.true_k_len;
                                            int gap = estimated_dist - bases_dist;

                                            if(bases_dist < estimated_dist + MAX_VARIANCE){
                                                contig_link_matrix[contig_a_indx][contig_b_indx].weight++;
                                                contig_link_matrix[contig_a_indx][contig_b_indx].total_dist += gap;
                                            }
                                        }
                                    }
                                    
                                }
                            }
                        }
                    }
                        
                }



                int find_intermediate_contig(
                    int target,
                    const vector<int>&possible_contigs,
                    int estimated_gap
                ){

                    //cout<<"finding middle contig\n";

                    if(estimated_gap < graph.true_k_len){return-1;}
                    int best_weight = MIN_WEIGHT;
                    int best_indx = -1;
                    vector<int>valid_links;
                    for(int c_i:possible_contigs){
                        if(contig_link_matrix[c_i][target].weight > best_weight){
                            if(contigs[c_i].nodes.size() * graph.true_k_len < estimated_gap + MAX_VARIANCE){
                                best_indx = c_i;
                                best_weight = contig_link_matrix[c_i][target].weight;
                            }
                        }
                    }
                    if(best_indx != -1){
                        int gap = contig_link_matrix[best_indx][target].dist();
                        if(gap > graph.true_k_len){
                            for(int j = 0; j<contig_link_matrix[best_indx].size(); j++){
                                if(
                                    contig_link_matrix[best_indx][j].weight >= MIN_WEIGHT
                                    &&
                                    contig_link_matrix[best_indx][j].dist() < graph.true_k_len
                                ){
                                    valid_links.push_back(j);
                                }
                            }
                            int intermediate = find_intermediate_contig(target,valid_links,gap);
                            if(intermediate != -1){
                                join_contigs(best_indx,intermediate);
                            }
                        }
                    }
                    return best_indx;
                }




                int find_best_contig_match(int i){
                    //cout<<"finding contig match\n";
                    int best_indx = -1;
                    int best_weight = MIN_WEIGHT;
                    int best_dist = 0;
                    vector<int>valid_links;
                    for(int j = 0; j<contig_link_matrix[i].size(); j++){
                        if(
                            contig_link_matrix[i][j].weight >= MIN_WEIGHT
                            &&
                            contig_link_matrix[i][j].dist() < graph.true_k_len
                        ){
                            valid_links.push_back(j);
                        }
                        if(contig_link_matrix[i][j].weight > best_weight){
                            best_indx = j;
                            best_weight = contig_link_matrix[i][j].weight;
                            best_dist = contig_link_matrix[i][j].dist();
                        }
                    }
                    if(best_indx != -1){
                        if(contig_link_matrix[i][best_indx].dist() > graph.true_k_len){
                            int middle_contig = find_intermediate_contig(best_indx,valid_links,best_dist);
                            if(middle_contig != -1){
                                join_contigs(i,middle_contig);
                            }
                        }
                    }
                    return best_indx;
                }


                void join_contigs(int i,int j){
                    cout<<"joined contigs "<< i<<" and "<<j<<"\n";
                    joined_contigs[j] = true;
                    contigs[i].insert(contigs[j]);
                }

                void clear_joined_contigs(){
                    for(int i = 0; i<contigs.size(); i++){
                        if(joined_contigs[i]){
                            contigs[i].nodes.clear();
                        }
                    }
                }


        };

};
    



    
int main(){
    vector<string>entries;
    string entry;
    vector<uint16_t>distances;
    bool is_read_pair = false;

    int count;
    cin>>count;
    for(int i = 0; i<count; i++){
        cin>>entry;
        int pos_1 = entry.find('|');
        if(pos_1 == string::npos){
            entries.push_back(entry);
        }
        else{
            is_read_pair = true;
            int pos_2 = entry.find('|',pos_1+1);
            string r1 = entry.substr(0,pos_1);
            string r2 = entry.substr(pos_1+1,pos_2-pos_1-1);
            entries.push_back(r1);
            entries.push_back(r2);
            uint16_t d = (std::stoi(entry.substr(pos_2+1)));
            distances.push_back(d);
        }
    }
    GENOME_ASSEMBLER assembler(entries,distances,is_read_pair);
    assembler.assemble_genome();
        
}