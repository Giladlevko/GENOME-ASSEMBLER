#pragma once

#ifndef DE_BRUIJN_GRAPH_HPP
#define DE_BRUIJN_GRAPH_HPP


#include "common.hpp"           //For STRING_REF
#include "kmer_bit_map.hpp"      //For K_MER_BIT_MAP

#include<vector>
#include<string>
#include<cstdint>
#include<queue>                 //For std::queue (used in TIP_REMOVER)
#include<unordered_set>         //For std::unordered_set
#include<unordered_map>         //For std::unordered_map (used in BUBBLE_REMOVER)
#include<stdexcept>             //For std::out_of_range


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



#endif