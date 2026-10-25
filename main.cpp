#include "source/genome_assembler.hpp"  //For GENOME_ASSEMBLER

#include<vector>
#include<string>
#include<cstdint>
#include<iostream>              //For std::cin, std::cout



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