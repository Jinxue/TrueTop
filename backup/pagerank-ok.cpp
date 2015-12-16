/**
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 
 *
 * @section DESCRIPTION
 *
 * Simple pagerank implementation. Uses the basic vertex-based API for
 * demonstration purposes. A faster implementation uses the functional API,
 * "pagerank_functional".
 */

#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>

#define GRAPHCHI_DISABLE_COMPRESSION


#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"
#include "graphchi_types.hpp"

using namespace graphchi;
 
//#define THRESHOLD 1e-7
#define RANDOMRESETPROB 0.0
//#define LOGOUTPUT
#define SEEDS

typedef float VertexDataType;
typedef float EdgeDataType;

struct PagerankProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {

    /*
     * We use two vectors to record the pagerank values in current and the last iterations 
     */
    std::vector<float> *last_it_pagerank;
    std::vector<float> *curr_it_pagerank;
    int num_vertices;
 
    // Define the top-K
    int K;

    // Error tolerance
    float tolerance;

    // termination type
    std::string terminate;

    // seed file
    std::map <vid_t, float> seedMap;

    PagerankProgram( int num_v, int k, float tol, std::string term, std::string seedF, std::string seedMode, unsigned int seedsN){
	num_vertices = num_v;
	std::cout << "The total verteces number is: " << num_vertices << std::endl;
	last_it_pagerank = new std::vector<float>(num_vertices, 0.0);
	curr_it_pagerank = new std::vector<float>(num_vertices, 0.0);
	K = k;
	tolerance = tol;
	terminate = term;
	buildSeeds(seedF, seedMode, seedsN);
    }    

    ~PagerankProgram(){
	last_it_pagerank->erase(last_it_pagerank->begin(), last_it_pagerank->end());
	delete last_it_pagerank;
	curr_it_pagerank->erase(curr_it_pagerank->begin(), curr_it_pagerank->end());
	delete curr_it_pagerank;
    }

    /**
     * Load the seed file with format: chiId ID rpg weight
     */
    void buildSeeds(std::string seedF, std::string mode, unsigned int seedsN){
	std::ifstream infile(seedF.c_str());
	vid_t chiid, id;
	float rpg, weight;
	float total_rpg = 0.0;

	unsigned int ii = 0;
	while(infile >> chiid >> id >> rpg >> weight){
	    if (ii < seedsN){
	        seedMap[chiid] = rpg;
	        total_rpg += rpg;
		ii ++;
	    } else
		seedMap[chiid] = 0.0;
	}

	if (seedsN > seedMap.size())
	    seedsN = seedMap.size();
	
	ii = 0;
	for (std::map<vid_t, float>::iterator it = seedMap.begin(); (it != seedMap.end() && (ii ++ <= seedsN)); it ++){
	    if (mode.compare("even") == 0)
		it->second = 1.0 / seedsN;
	    else if (mode.compare("rpg") == 0)
		it->second = it->second / total_rpg;
	    //std::cout << it->first << " : " << it->second << std::endl;
	}
	std::cout << " #Seeds : " << seedsN << std::endl;
    }

    /**
      * Called before an iteration starts. Not implemented.
      */
    void before_iteration(int iteration, graphchi_context &info) {
    }
    

    /*
     * Sort the vector with index
     */
    template <typename T>
    std::vector<int> sort_indexes(const std::vector<T> &v){
	// Initialize original index locations
	std::vector<int> idx(v.size(), 0);
	//std::cout << "list size: " << v.size() << std::endl;
	for (size_t i = 0; i < idx.size(); i ++)
	    idx[i] = i;
	
	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(), [&v](int i1, int i2) {
	    return v[i1] > v[i2];}
	);


	return idx;
    }

    /**
     * Obtain the rank for the sorted vector
     */
    template <typename T>
    std::vector<int> rank_indexes(const std::vector<T> &v){
	// obtain the index for the sorting
	std::vector<int> idx = sort_indexes(v);
	
	std::vector<int> rank(v.size(), 0);
	for (size_t i = 0; i < rank.size(); i ++)
	    rank[idx[i]] = i;
	    //rank[i] = idx[i];
	return rank;
    }

    /**
     * Get the error for the top K nodes
     */
    float topk_diff (){
	// Obtain the topK nodes for the last_it_pagerank
	std::vector<int> idx_last = sort_indexes(*last_it_pagerank);

	// For each of top-K nodes in the last_it_pagerank, find its ranking in the curr_it_pagerank
	std::vector<int> rank_curr = rank_indexes(*curr_it_pagerank);

	int err_sum = 0;
	for (int a = 0; a < K; a ++){
	    err_sum += std::abs(a - rank_curr[idx_last[a]]);
	    //std::cout << "error: " << err_sum << "rank : " << a  << "---"  <<rank_curr[idx_last[a]] << std::endl;
	}
	return err_sum;
    }

    /**
      * the difference between two iterations
      */
    float pagerank_diff() {
    	// We first check the error sum of the iteration
	float err_sum = 0;
	for (int i = 0; i < num_vertices; i ++){
#ifdef LOGOUTPUT
	    if (i == 1341)
    	        std::cout <<  i << " -> current: " << curr_it_pagerank->at(i) << ", last: " << last_it_pagerank->at(i) << std::endl;
#endif
	    err_sum += std::abs(last_it_pagerank->at(i) - curr_it_pagerank->at(i));
	    //last_it_pagerank->at(i) = curr_it_pagerank->at(i);
	}
	//return err_sum / num_vertices;
	return err_sum;
    }
    
    /**
      * Called after an iteration has finished. Not implemented.
      */
    void after_iteration(int iteration, graphchi_context &ginfo) {
	
	if (terminate.compare("topK") == 0){
            float err_sum = 0;
	
	    // Terminate it by the ranking of top K
	    err_sum = topk_diff();
	    std::cout << "At the iteration " << iteration << " with error: " << err_sum << std::endl;
	    if ((iteration > 0) && (err_sum <= tolerance)) {
	        std::cout << "Achieve the error tolerace, terminate." << std::endl;
	        ginfo.set_last_iteration(ginfo.iteration);
	    }
	} else if (terminate.compare("total") == 0){
	    // Terminate it by the error tolerance
	    float err_sum = pagerank_diff();
//#ifdef LOGOUTPUT
	    std::cout << "At the iteration " << iteration << " with error: " << err_sum << std::endl;
//#endif
	    if (err_sum < tolerance){
	        std::cout << "Achieve the error tolerace, terminate." << std::endl;
	        ginfo.set_last_iteration(ginfo.iteration);
	    }
	} else{
	    std::cout<< "Termination doesn't support: " << terminate << std::endl;
	    return;
	}

	/* Finally, we update the last_it_pagerank */
	for (int i = 0; i < num_vertices; i ++){
	    last_it_pagerank->at(i) = curr_it_pagerank->at(i);
	}
	
    }

    /**
      * Called before an execution interval is started. Not implemented.
      */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
    }
    
    
    /**
      * One-time update the weight of the edge.
      * The weight of each edge is the normalized ratio from the source node
      */
    void update_weight(graphchi_vertex<VertexDataType, EdgeDataType> &v){
	float sum = 0.0;
        for(int i=0; i < v.num_outedges(); i++) {
            graphchi_edge<EdgeDataType> * edge = v.outedge(i);
	    //We sum the weight values for each outgoing edges.
	    float weight = edge->get_data();
	    sum += weight;
	    //sum ++;
        }
	
        for(int i=0; i < v.num_outedges(); i++) {
            graphchi_edge<EdgeDataType> * edge = v.outedge(i);
	    float weight = edge->get_data();
	    float norm_weight = 1.0 * weight / sum;
	    //eData.pagerank = 1.0 * quota / sum;
	    edge->set_data(norm_weight);
#ifdef LOGOUTPUT
	    if (v.id() == 1341 || edge->vertex_id() == 1341)
    	        std::cout <<  v.id() << " -> " << edge->vertex_id() << " with data: " << eData.pagerank << " with weight: " << eData.weight << "quota:" << quota << std::endl;
#endif
         }
    }


    /**
      * Pagerank update function.
      */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {
        float sum = 0.0;
	float pagerank = 0.0;
        if (ginfo.iteration == 0) {
            /* On first iteration, initialize vertex and out-edges. 
               The initialization is important,
               because on every run, GraphChi will modify the data in the edges on disk. 
             */
#ifdef SEEDS
	    // For the seeds
	    if((v.id() > 0) && (seedMap.find(v.id()) != seedMap.end()))
            	pagerank = seedMap.find(v.id()) -> second;
	    else
		pagerank = 0.0;
#else
	    pagerank = 1.0 / (num_vertices - 1);
#endif
	    //std::cout << v.id() << "'s pagerank: " << pagerank << std::endl;
	    // Update the weigh for the edge
	    update_weight(v);

        } else {
            /* Compute the sum of neighbors' weighted pageranks by
               reading from the in-edges. */
            for(int i=0; i < v.num_inedges(); i++) {
		graphchi_edge<EdgeDataType> * edge = v.inedge(i);
		sum += (last_it_pagerank->at(edge->vertex_id())) * (edge->get_data());
#ifdef LOGOUTPUT
	        if (v.id() == 1341)
    	            std::cout <<  "Update: " << v.inedge(i)->vertex_id() << "->" << v.id() << " sum: " << sum << std::endl;
#endif
            }
            
            /* Compute my pagerank */
            pagerank = RANDOMRESETPROB + (1.0 - RANDOMRESETPROB) * sum;
        }

        /* Set my new pagerank as the vertex value */
        v.set_data(pagerank); 
	
	/* Update the current pagerank vector */
	curr_it_pagerank->at(v.id()) = pagerank;
#ifdef LOGOUTPUT
	if (v.id() == 1341)
    	    std::cout <<  v.id() << " -> " << curr_it_pagerank->at(v.id()) << std::endl;
#endif
    }
    
};

int main(int argc, const char ** argv) {
    graphchi_init(argc, argv);
    metrics m("pagerank");
    
    /* Parameters */
    std::string filename    = get_option_string("file"); // Base filename
    int niters              = get_option_int("niters", 4);
    bool scheduler          = false;                    // Non-dynamic version of pagerank.
    int ntop                = get_option_int("top", 20);
    int topk                = get_option_int("topk", 100);
    float tolerance         = get_option_float("tolerance", 100);
    std::string terminate   = get_option_string("terminate");
    std::string seedF 	    = get_option_string("seedsfile", "seeds.txt");
    std::string seedMode    = get_option_string("seedsmode", "even");
    unsigned int seedsN    = get_option_int("seedsnum", 100);
    
    /* Process input file - if not already preprocessed */
    int nshards             = convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));

    /* Run */
    //graphchi_engine<float, float> engine(filename, nshards, scheduler, m); 
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
    engine.set_modifies_inedges(false); // Improves I/O performance.
    std::cout << get_num_vertices(filename) << " " << topk << " " << tolerance << " " << terminate << " " << seedF << " " << seedMode << " " << seedsN << std::endl;
    PagerankProgram program(get_num_vertices(filename), topk, tolerance, terminate, seedF, seedMode, seedsN);
    engine.run(program, niters);
        
    /* Output top ranked vertices to a file */
    std::ofstream ofile;
    std::ostringstream outfName;
    //outfName << filename << "-pagerank-" << topk << "-" << tolerance << ".txt";
    // For the seed evaluation
    outfName << filename << "-pagerank-seeds-" << seedMode << "-" << seedsN << ".txt";
    //outfName << filename << "-pagerank" << ".txt";
    ofile.open(outfName.str().c_str());
    std::vector< vertex_value<float> > top = get_top_vertices<float>(filename, ntop);
    std::cout << "Print top " << ntop << " vertices:" << std::endl;
    //size_t num_vertices = get_num_vertices(filename);
    for(int i=0; i < (int)top.size(); i++) {
        if(i < 30)
            std::cout << (i+1) << ". " << top[i].vertex << "\t" << top[i].value << std::endl;
        ofile << top[i].vertex << "\t" << top[i].value << std::endl;
    }
    
    metrics_report(m);    
    return 0;
}

