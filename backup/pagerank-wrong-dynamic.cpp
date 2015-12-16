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
#include <iostream>
#include <fstream>
#include <cmath>

#define GRAPHCHI_DISABLE_COMPRESSION
#define DYNAMICEDATA 1


#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"
#include "api/dynamicdata/chivector.hpp"
#include "api/chifilenames.hpp"

using namespace graphchi;
 
//#define THRESHOLD 1e-1    
#define RANDOMRESETPROB 0.0


typedef float VertexDataType;
#ifndef DYNAMICEDATA
typedef float EdgeDataType;
#else
typedef chivector<float>  EdgeDataType;
#endif

struct PagerankProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    /**
      * Called before an iteration starts. Not implemented.
      */
    void before_iteration(int iteration, graphchi_context &info) {
    }
    
    /**
      * Called after an iteration has finished. Not implemented.
      */
    void after_iteration(int iteration, graphchi_context &ginfo) {
    }
    
    /**
      * Called before an execution interval is started. Not implemented.
      */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
    }
    
    /**
      * Update the weigthed edge chivector
      * We first obtain the edge weight from the first element, sum them, then update the 
      * second item by eacg edge's weight
      */
    void update_edge_data(graphchi_vertex<VertexDataType, EdgeDataType> &v, float quota, bool first){
	float sum = 0.0;
	//if(first)
            for(int i=0; i < v.num_outedges(); i++) {
                graphchi_edge<EdgeDataType> * edge = v.outedge(i);
                if (edge != NULL) {
                    chivector<float> * evector = edge->get_vector();
		    //std::cout << evector->size() << std::endl;
		    /*if (first)
                        assert(evector->size() == 1);
		    else
                        assert(evector->size() == 2);
                    assert(evector->size() == 2);*/
    	            std::cout <<  v.id() << " with data: " << evector->get(0) << std::endl;
	            sum += evector->get(0);
		    /*if (first){
                        evector->add(sum);
                	assert(evector->size() == 2);
		    }*/
	        }
            }
	
        for(int i=0; i < v.num_outedges(); i++) {
            graphchi_edge<EdgeDataType> * edge = v.outedge(i);
            if (edge != NULL) {
                chivector<float> * evector = edge->get_vector();
//                assert(evector->size() == 2);
		float val = quota * evector->get(0) / sum;
		//evector->set(1, val);
		if(first && (evector->size() == 1))
                        evector->add(val);	
		evector->set(1, val);
    		//std::cout <<  v.id() << " with data: " << evector->get(0) << std::endl;
	    }
         }
    }

    /**
      * Pagerank update function.
      */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {
        float sum=0;
	float prv = 0.0;
	float pagerankcont = 0.0;	

        if (ginfo.iteration == 0) {
            /* On first iteration, initialize vertex and out-edges. 
               The initialization is important,
               because on every run, GraphChi will modify the data in the edges on disk. 
             */
	    /* For the weighted version */
	    update_edge_data(v, 1.0, true);
            v.set_data(RANDOMRESETPROB); 
            //v.set_data(1.0); 
        } else {
	    /* We need to come up with the weighted version */
            for(int i=0; i < v.num_inedges(); i++) {
                chivector<float> * evector = v.inedge(i)->get_vector();
                assert(evector->size() >= 2);
                sum += evector->get(1);
    		//std::cout <<  v.id() << " with data: " << evector->get(1) << " with weight " << evector->get(0) << std::endl;
    		//std::cout <<  v.id() << " edge endpoint: " << v.inedge(i)->vertex_id() << std::endl;
		//evector->clear();
	    }

            /* Compute my pagerank */
            prv = RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum;
	    //std::cout << "sum" << sum << "pagerank: " << prv << std::endl;

	    update_edge_data(v, prv, false);
            /* Keep track of the progression of the computation.
               GraphChi engine writes a file filename.deltalog. */
	    double delta = std::abs(prv - v.get_data());
	    //std::cout << "pagerank: " << prv << "v.data" << v.get_data() << "delta: " << delta << std::endl;
            ginfo.log_change(delta);
            
            /* Set my new pagerank as the vertex value */
            v.set_data(prv);
        }
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
    
    /* Process input file - if not already preprocessed */
    int nshards             = convert_if_notexists<float>(filename, get_option_string("nshards", "auto"));

    /* Run */
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m);
    engine.set_modifies_inedges(false); // Improves I/O performance.
    PagerankProgram program;
    engine.run(program, niters);
        
    /* Output top ranked vertices */
    std::ofstream ofile;
    ofile.open((filename + "-pagerank.txt").c_str());
    std::vector< vertex_value<float> > top = get_top_vertices<float>(filename, ntop);
    std::cout << "Print top " << ntop << " vertices:" << std::endl;
    size_t num_vertices = get_num_vertices(filename);
    for(int i=0; i < (int)top.size(); i++) {
        if(i < 30)
            std::cout << (i+1) << ". " << top[i].vertex << "\t" << top[i].value / num_vertices << std::endl;
        ofile << top[i].vertex << "\t" << top[i].value /num_vertices << std::endl;
    }
    
    metrics_report(m);    
    ofile.close();
    return 0;
}

