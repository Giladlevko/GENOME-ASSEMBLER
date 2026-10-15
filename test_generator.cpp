#include<vector>
#include<random>
#include<string>
#include<fstream>
#include<algorithm>
#include<iostream>

using namespace std;

const bool genarate_from_real_gene = true;

enum class GENE_TYPES{N_DELTO,E_COLI};

const GENE_TYPES GENE_TYPE = GENE_TYPES::N_DELTO;

const int gene_len = 4600000;
const int GC_percentage = 50;

const int read_len_fluctuation = 50;
const int read_len = 100;

const int max_read_len = read_len + read_len_fluctuation;
const int min_read_len = read_len - read_len_fluctuation;

const int coverage = 30;

const string GC_bases = "CG";
const string AT_bases = "AT";
const string bases = "ACGT";

const bool has_errors = true;

const bool is_read_pair = true;
const int pair_dist_max = 300;
const int pair_dist_min = 100;
const int dist_variance = 20;


char get_rand_base(mt19937& rng){
    uniform_int_distribution<int>rand_perc(0,100);
    uniform_int_distribution<int>rand_index(0,1);
    if(rand_perc(rng)<=GC_percentage){
        return GC_bases[rand_index(rng)];
    }
    else{
        return AT_bases[rand_index(rng)];
    }

}


void gen_rand_genome(string& gene){
    
    mt19937 rng(100);
    uniform_int_distribution<int>dist(0,3);
    for(int i = 0; i<gene_len; i++){
        gene += get_rand_base(rng);
    }
}

string get_gene_file(){
    switch(GENE_TYPE){
        case GENE_TYPES::N_DELTO: return("genome_assembly/N_deltocephalinicola_gene.txt");
        case GENE_TYPES::E_COLI: return("genome_assembly/ecoli_gene.txt");
        default: return("genome_assembly/ecoli_gene.txt");
    }
}


void read_real_gene(string&gene){
    string path = get_gene_file();
    ifstream file(path);
    string read;
    string name;
    file>>name;
    cout<<"GENE NAME: "<<name<<"\n";
    while(file>>read){
        gene+=read;
    }
}


void diff_rand_base(char&base,mt19937& rng){
    string shuffled_bases = bases;
    shuffle(shuffled_bases.begin(),shuffled_bases.end(),rng);
    for(const char& c: shuffled_bases){
        if(c != base){
            base = c;
            return;
        }
    }
}

void break_gene_into_reads(const string& gene,vector<string>&reads){
    mt19937 rng(100);

    uniform_int_distribution<int>read_fluctuation(
        min_read_len,
        max_read_len
    );

    
    int read_count = gene.size() * coverage / read_len;
    if(is_read_pair){
        read_count /= 2;
    }

    for(int i = 0; i<read_count; i++){

        int actual_d;
        int mean_d;
        
        int actual_read_len = read_fluctuation(rng);

        int max_start = gene.size()-actual_read_len;
        int actual_max_start = max_start;

        if(is_read_pair){
            uniform_int_distribution<int>mean_d_dist(pair_dist_min,pair_dist_max);
            mean_d = mean_d_dist(rng);
            uniform_int_distribution<int>d_dist(mean_d - dist_variance, mean_d + dist_variance);
            actual_d = d_dist(rng);
            actual_max_start = max_start - (actual_d + actual_read_len);
        }

        uniform_int_distribution<int>dist(0,actual_max_start);

        int start = dist(rng);
        string read = gene.substr(start,actual_read_len);
        uniform_int_distribution<int>read_dist(0,read.size()-1);
        if(has_errors){
            //introduce a 1% error
            int error_i = read_dist(rng);
            diff_rand_base(read.at(error_i),rng);
        }
        

        if(is_read_pair){
            int r2_start = start + actual_read_len + actual_d;
            string read_2 = gene.substr(r2_start, actual_read_len);
            if(has_errors){
                //introduce a 1% error
                int error_i = read_dist(rng);
                diff_rand_base(read_2.at(error_i),rng);
            }

            read += '|' + read_2 + '|' + to_string(mean_d + actual_read_len);
        }
        reads.push_back(read);

    }
}

struct my_bool{
    my_bool(bool t):is_true(t){}
    bool is_true;
};

ostream& operator<<(ostream&out,my_bool is_true){
    return out << (is_true.is_true ? "True" : "False"); 
}

float actual_gc(const string&gene){
    double gc_count = 0;
    for(const char& base:gene){
        if(base == 'C' || base == 'G'){
            gc_count+=1.0;
        }
    }
    return (gc_count / gene.size()) * 100.0;
}

void write_test_to_file(){
    string gene = "";
    if(!genarate_from_real_gene){
        gen_rand_genome(gene);
    }
    else{
        read_real_gene(gene);
    }
    ofstream gene_file("genome_assembly/reference_gene.txt");
    gene_file<<">REFERENCE GENOME\n"<<gene;
    vector<string>reads;
    break_gene_into_reads(gene,reads);
    ofstream file("genome_assembly/test_inputs.txt");
    if(!file){cout<<"could not open file\n";}

    cout<<"Writing to file: genome_assembly/test_inputs.txt\nGenome size: "<<gene.size()
    <<"\nAverage read size: "<<read_len<<"\nRead size fluctuation: "<<min_read_len<<" - "<<max_read_len
    <<"\nRead count: "<<reads.size()<<"\nError percentage: 1%"
    <<"\nCoverage: "<<coverage <<"\nPaired inputs: "<<my_bool(is_read_pair)<<"\nGC% "<< actual_gc(gene)
    <<"\nHas errors: "<<my_bool(has_errors);
    
    file<<reads.size()<<"\n";
    
    for(int i = 0; i<reads.size(); i++){
        file<<reads[i]<<"\n";
    }
    file.close();
}

int main(){
    write_test_to_file();
    return 0;
}