#pragma once

#ifndef KMER_MAP_HPP
#define KMER_MAP_HPP

#include "common.hpp"       //For STRING_REF

#include<vector>
#include<string>
#include<unordered_map>
#include<algorithm>         //For std::lower_bound, std::min, std::max, std::sort, std::distance
#include<unordered_set>
#include<utility>           //For std::pair
#include<stdexcept>         //For std::out_of_range
#include<cstdint>
#include<iostream>          //For std::ostream, std::cout




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


#endif