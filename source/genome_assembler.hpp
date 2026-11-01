#pragma once

#ifndef GENOME_ASSEMBLER_HPP
#define GENOME_ASSEMBLER_HPP

#include "de_bruijn_graph.hpp"  // For DE_BRUIJN_GRAPH, K_MER_BIT_MAP, and STRING_REF

#include <vector>
#include <string>
#include <cstdint>
#include <fstream>              // For std::ofstream (used in assemble_genome)
#include <iostream>             // For std::cout, std::endl, std::ostream
#include <stdexcept>            // For std::out_of_range



    
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
            cout<<"\nGenome assembled!\ncontigs outputed at "<<contig_file_path<<"\a"<<endl;
            
        }

        
    private:
        
        DE_BRUIJN_GRAPH graph;

        bool is_read_pair;
        vector<uint16_t>&distances;
        vector<string>&reads;


        struct CONTIG{
            CONTIG(vector<int>a):nodes(a){}
            vector<int>nodes;
            void insert(const CONTIG&other,int gap){
                //cout<<"gap is "<<gap<<" main size = "<<size()<<" other size = "<<other.size()<<"\n";
                if(gap <= 0){
                    nodes.insert(nodes.end(),other.nodes.begin()+abs(gap),other.nodes.end());
                }
                else{
                    vector<int>n_bases(gap,-1);
                    nodes.insert(nodes.end(),n_bases.begin(),n_bases.end());
                    nodes.insert(nodes.end(),other.nodes.begin(),other.nodes.end());
                }
               
            }
            int& operator[](int i){
                if(i > nodes.size()){
                    throw std::out_of_range("contig base index out of range");
                }
                return nodes[i];
            }
            size_t size()const{
                return (nodes.size());
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

            
            const int min_contig_len = 0;

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
                    if(contigs[i][j] == -1){
                        out<<"N";
                        //cout<<"new unknown base\n";
                    }
                    else{
                        id_to_str.print_back_char_at(contigs[i][j],out);
                    }
                    
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
                    
                }

                void join_linked_contigs(){
                    for(int i = 0; i<5; i++){
                        assign_node_to_contig_index();
                        build_contig_link_matrix();
                        //node to contig is not needed anymore so I can clear it
                        node_to_contig.clear(); node_to_contig.shrink_to_fit();
                        for(int i = 0; i<contig_link_matrix.size(); i++){
                            if(joined_contigs[i]){continue;}
                            int estimated_gap = 0;
                            int c = find_best_contig_match(i,estimated_gap);
                            if(c != -1 && contig_link_is_mutual_best(i,c)){
                                join_contigs(i,c,estimated_gap);
                            }
                        
                        }
                        MIN_WEIGHT += 5;
                        //MAX_VARIANCE += 1;
                        clear_joined_contigs();
                    }
                    
                    cout<<"finished linking\n";
                    
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

                int MIN_WEIGHT = 3;
                int MAX_VARIANCE = 30;
                int MAX_GAP = 120;

                void assign_node_to_contig_index(){
                    node_to_contig.assign(graph.size(),{});
                    for(int i = 0; i<contigs.size(); i++){
                        for(int j = 0; j<contigs[i].nodes.size(); j++){
                            int v = contigs[i].nodes[j];
                            if(v != -1){
                                node_to_contig[v].push_back({i,j});
                            }
                            
                        }
                    }
                }


                void build_contig_link_matrix(){
                    cout<<"building links\n";
                    contig_link_matrix.assign(contigs.size(),vector<CONTIG_EDGE>(contigs.size()));
                    for(auto it:graph.id_to_str.pair_map){
                        if(node_to_contig[it.first].size() > 1){continue;}
                        for(int i = 0; i<node_to_contig[it.first].size(); i++){
                            int contig_a_indx = node_to_contig[it.first][i].contig_index;
                            if(contig_a_indx>=0){

                                int contig_a_size = contigs[contig_a_indx].nodes.size();
                                int a_node_depth = node_to_contig[it.first][i].depth;

                                int a_dist_from_end =  contig_a_size - a_node_depth; 

                                for(const auto& data:it.second){
                                    if(node_to_contig[data.target_index].size() > 1){continue;}
                                    for(int j = 0; j<node_to_contig[data.target_index].size(); j++){
                                        int contig_b_indx = node_to_contig[data.target_index][j].contig_index;
                                        if(contig_b_indx >= 0){
                                            int b_node_depth = node_to_contig[data.target_index][j].depth;
                                            int estimated_dist = data.dist;
                                            int bases_dist = (a_dist_from_end + b_node_depth) + graph.true_k_len;
                                            int gap = estimated_dist - bases_dist;

                                            //TO DO currently this only blocks large
                                            //gaps but should also block large overlaps
                                            //simply doing abs causes missasemblies for some reason
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
                    int estimated_gap,
                    int search_depth = 0
                ){

                    //cout<<"finding middle contig\n";

                    if(estimated_gap < MAX_GAP || search_depth > 3){return-1;}
                    int min_middle_weight = min(1,MIN_WEIGHT/3);
                    int best_weight = min_middle_weight;
                    int best_indx = -1;
                    vector<int>valid_links;
                    for(int c_i:possible_contigs){
                        if(contig_link_matrix[c_i][target].weight > best_weight){
                            if(contigs[c_i].nodes.size() + graph.true_k_len-1 < estimated_gap + MAX_VARIANCE){
                                best_indx = c_i;
                                best_weight = contig_link_matrix[c_i][target].weight;
                            }
                        }
                    }
                    if(best_indx != -1){
                        int gap = contig_link_matrix[best_indx][target].dist();
                        if(gap > MAX_GAP){
                            for(int j = 0; j<contig_link_matrix[best_indx].size(); j++){
                                if(
                                    contig_link_matrix[best_indx][j].weight >= min_middle_weight
                                    &&
                                    contig_link_matrix[best_indx][j].dist() < MAX_GAP
                                ){
                                    valid_links.push_back(j);
                                }
                            }
                            int intermediate = find_intermediate_contig(target,valid_links,gap,++search_depth);
                            if(intermediate != -1){
                                int intermediate_gap = contig_link_matrix[best_indx][intermediate].dist();
                                join_contigs(best_indx,intermediate,intermediate_gap);
                            }
                            else{
                                best_indx = -1;
                            }
                        }
                    }
                    return best_indx;
                }




                int find_best_contig_match(int i,int& gap){
                    //cout<<"finding contig match\n";
                    int best_indx = -1;
                    int best_weight = MIN_WEIGHT;
                    int second_weight = 0;
                    int best_dist = 0;
                    vector<int>valid_links;
                    for(int j = 0; j<contig_link_matrix[i].size(); j++){
                        if(
                            contig_link_matrix[i][j].weight >= MIN_WEIGHT
                            &&
                            contig_link_matrix[i][j].dist() < MAX_GAP
                        ){
                            valid_links.push_back(j);
                        }
                        if(contig_link_matrix[i][j].weight > best_weight){
                            best_indx = j;
                            second_weight = best_weight;
                            best_weight = contig_link_matrix[i][j].weight;
                            best_dist = contig_link_matrix[i][j].dist();
                        }
                        else if(second_weight < contig_link_matrix[i][j].weight){
                            second_weight = contig_link_matrix[i][j].weight;
                        }
                    }
                    if(best_indx != -1){
                        if(contig_link_matrix[i][best_indx].dist() > MAX_GAP){
                            int middle_contig = find_intermediate_contig(best_indx,valid_links,best_dist);
                            if(middle_contig != -1){
                                int middle_gap = contig_link_matrix[i][middle_contig].dist();
                                cout<<"middle gap "<<middle_gap<<"\n";
                                join_contigs(i,middle_contig,middle_gap);
                                cout<<"gap changed from "<<best_dist;
                                best_dist -= contigs[middle_contig].size();// + graph.true_k_len - 1;
                                cout<<" to "<<best_dist<<"\n";
                            }
                            else{
                                best_indx = -1;
                                
                            }
                        }
                    }
                    //if(best_weight < second_weight * 1.8){best_indx = -1;}
                    gap = best_dist;
                    return best_indx;
                }

                //TO DO add a second best and choose the link only if 
                //the best is much higher than the second best
                //I will have to test how much higher
                bool contig_link_is_mutual_best(int i, int j){
                    int best_pre = -1;
                    int best_weight = 0;
                    int second_weight = 0;
                    for(int k = 0; k<contig_link_matrix.size(); k++){
                        if(contig_link_matrix[k][j].weight>best_weight){
                            second_weight = best_weight;
                            best_weight = contig_link_matrix[k][j].weight;
                            best_pre = k;
                        }
                        else if(second_weight < contig_link_matrix[k][j].weight){
                            second_weight = contig_link_matrix[k][j].weight;
                        }
                    }
                    /*
                    return ( 
                        contig_link_matrix[i][j].weight >= contig_link_matrix[best_pre][j].weight
                        &&
                        best_weight > second_weight * 1.5 
                    );
                    */
                    //return ( contig_link_matrix[i][j].weight >= contig_link_matrix[best_pre][j].weight*0.9 );
                    return ( i == best_pre );
                }


                void join_contigs(int i,int j,int gap = 0){
                    //cout<<"joined contigs "<< i<<" and "<<j<<"\n";
                    if(gap<0 && abs(gap)>= contigs[j].size()-1){return;}
                    cout<<"gap is "<<gap<<"\n";
                    joined_contigs[j] = true;
                    contigs[i].insert(contigs[j],gap);
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


#endif
