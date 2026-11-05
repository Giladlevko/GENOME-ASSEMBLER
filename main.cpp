#include "source/genome_assembler.hpp"  //For GENOME_ASSEMBLER

#include<vector>
#include<string>
#include<cstdint>
#include<fstream>
#include<iostream>              //For std::cin, std::cout



int main(int argc, char* argv[]){
    string contig_out_path = "genome_assembly/contig_output.txt";
    string test_input_path = "genome_assembly/test_inputs.txt";
    if(argc < 3){
        cout<<"Not enough files inputed!\n" 
        <<"Usage: ./main <input_file> <output_file>\n"
        <<"Using default files!\n"
        <<"Input: genome_assembly/test_inputs.txt\n"
        <<"Contig output: genome_assembly/contig_output.txt\n";
    }
    else{
        test_input_path = argv[1];
        contig_out_path = argv[2];
    }

    vector<string>entries;
    string entry;
    vector<uint16_t>distances;
    bool is_read_pair = false;

    ifstream in_file(test_input_path);

    int count;
    in_file>>count;
    for(int i = 0; i<count; i++){
        in_file>>entry;
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
    GENOME_ASSEMBLER assembler(entries,distances,is_read_pair,contig_out_path);
    assembler.assemble_genome();
        
}