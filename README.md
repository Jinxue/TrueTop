This project is based on GraphChi platform from CMU (Thanks very much!).

My major job is to add an application function https://github.com/Jinxue/TrueTop/blob/master/example_apps/pagerank.cpp
to early terminate the PageRank by monitoring the top-K stablity, and some supporting files.

For the algorithm details, please check my TON paper: http://arxiv.org/abs/1504.04663,.

####Run-time performance:
On our desktop with 3.4GHz Intel-i7 3770 CPU, 16G Memory, a 7200RPM hard disk, and Ubuntu 12.04 LTS, one single iteration of credit distribution took 0.3s, 2.5s, 9.2s, and 17.1s for our four datasets in Table I with 4K, 10K, 1M and 2M nodes, respectively. For a graph with 2M nodes, TrueTop can thus find the top-1000 influential users after 1,000 iterations within less than five hours on
a commodity PC. Since TrueTop is expected to be run by a service provider with much more powerful computation resources, its runtime performance should be acceptable.


# GraphChi - disk-based large-scale graph computation

**NOTE: This project has been recently moved from Google Code, and some of the wiki pages might be partly broken.**

MIT Technology Review article about GraphChi: ["Your laptop can now analyze big data"](http://www.technologyreview.com/news/428497/your-laptop-can-now-analyze-big-data/?nlid=nldly&nld=2012-07-17)

### NEW: Graph contraction algorithms

Read about graph contraction technique, which we used to implemented efficient minimum spanning forest computation [ Graph Contraction Algorithms ](https://github.com/GraphChi/graphchi-cpp/wiki/Graph-Contraction-Algorithms).

### NEW: GraphChi version 0.2

Highlight of the release includes support for dynamically sized edge and vertex values.

**Read [release announcement.](https://github.com/GraphChi/graphchi-cpp/wiki/GraphChi-Version0p2-Release)  **

### Discussion group

http://groups.google.com/group/graphchi-discuss

## Introduction

GraphChi is a spin-off of the GraphLab ( http://www.graphlab.org ) -project from the Carnegie Mellon University. It is based on research by Aapo Kyrola ( http://www.cs.cmu.edu/~akyrola/) and his advisors. 

GraphChi can run very large graph computations on just a single machine, by using a novel algorithm for processing the graph from disk (SSD or hard drive). Programs for GraphChi are written in the *vertex-centric* model, proposed by GraphLab and Google's Pregel. GraphChi runs vertex-centric programs **asynchronously** (i.e changes written to edges are immediately visible to subsequent computation), and in **parallel**.  GraphChi also supports **streaming graph updates** and removal of edges from the graph. Section 'Performance' contains some examples of applications implemented for GraphChi and their running times on GraphChi.

The promise of GraphChi is to bring web-scale graph computation, such as analysis of social networks, available to anyone with a modern laptop. It saves you from the hassle and costs of working with a distributed cluster or cloud services. We find it much easier to debug applications on a single computer than trying to understand how a distributed algorithm is executed. 

In some cases GraphChi can solve bigger problems in reasonable time than many other available *distributed* frameworks. GraphChi also runs efficiently on servers with plenty of memory, and can use multiple disks in parallel by striping the data.

Even if you do require the processing power of high-performance clusters, GraphChi can be an excellent tool for developing and debugging your algorithms prior to deploying them to the cluster. For high-performance graph computation in the distributed setting, we direct you to GraphLab's new version (v2.1), which can now handle large graphs in astonishing speed. GraphChi supports also most of the new GraphLab v2.1 API (with some restrictions), making the transition easy.

GraphChi is implemented in plain C++, and available as open-source under the flexible Apache License 2.0.

### Java version

Java-version of GraphChi: https://github.com/GraphChi/graphchi-java

### Publication

GraphChi is part of the OSDI'12 proceedings. PDF of the paper can be downloaded here: http://select.cs.cmu.edu/publications/paperdir/osdi2012-kyrola-blelloch-guestrin.pdf

Slides (OSDI talk): http://www.cs.cmu.edu/~akyrola/files/osditalk-graphchi.pptx

### Collaborative Filtering Toolkit

Danny Bickson has ported several collaborative filtering algorithms from GraphLab to GraphChi: 
http://bickson.blogspot.com/2012/12/collaborative-filtering-with-graphchi.html

## Features

- Vertex-centric computation model (similar to GraphLab, Pregel or Giraph)
  ** Wrapper for GraphLab 2.1 API (Gather-Apply-Scatter model)
- Asynchronous, parallel execution, with (optional) deterministic scheduling (see semantics section at [Creating-GraphChi-Applications](https://github.com/GraphChi/graphchi-cpp/wiki/Creating-GraphChi-Applications) )

- Can run graphs with billions of edges, with linear scalability, on a standard consumer grade machine
    ** Can also utilize large amounts of memory by preloading (caching), making it competitive on large servers: UsingMoreMemory
- Multidisk striping - RAID-style operation, see MultipleDisksSupports
- Works well on both hard-drive and SSD.
  
- Evolving graphs, streaming graph updates (can add and delete edges): https://github.com/GraphChi/graphchi-cpp/wiki/Evolving-And-StreamingGraphs 

- Easy to install, headers-only, no dependencies.
- Tested on Mac OS X and Linux

## Getting Started

Best way to get started is to start from the ExampleApps page.
Prior to that, you need to download the source code (no configuration
or installation is required).

For an introduction on writing your own applications, read  [Creating-GraphChi-Applications](https://github.com/GraphChi/graphchi-cpp/wiki/Creating-GraphChi-Applications).

## How GraphChi works

GraphChi is based on the Parallel Sliding Windows method which allows efficient asynchronous processing of mutable graphs from disk. See [Introduction-To-GraphChi](https://github.com/GraphChi/graphchi-cpp/wiki/Introduction-To-GraphChi) for description.

## Performance

In the table below, we have picked some recent running time results for large-scale graph problems from the literature, and run the same experiment using GraphChi. For GraphChi, we used a Mac Mini (2012 model), with 8 gigabytes of RAM and 256 gigabyte SSD drive. 

While distributed clusters can solve the same problems faster than GraphChi on a single computer, for many purposes GraphChi's performance should be adequate. The numbers below do not include time for transferring the input to cloud or cluster, and usually do not include the graph loading time. GraphChi's running times include loading the graph and saving the results, but not preprocessing time. Preprocessing needs to be done only once per graph (you can run many different algorithms on the same preprocessed graph). The preprocessing times are listed in a separate table below.


<table class="wikitable"><tr><td style="border: 1px solid #ccc; padding: 5px;"> <strong>Application</strong> </td><td style="border: 1px solid #ccc; padding: 5px;"> <strong>Input graph</strong> </td><td style="border: 1px solid #ccc; padding: 5px;"> <strong>Graph size</strong></td><td style="border: 1px solid #ccc; padding: 5px;"> <strong>Comparison</strong> </td><td style="border: 1px solid #ccc; padding: 5px;"> <strong>GraphChi on Mac Mini (SSD)</strong> </td><td style="border: 1px solid #ccc; padding: 5px;"> Ref </td></tr> <tr><td style="border: 1px solid #ccc; padding: 5px;"> Pagerank - 3 iterations </td><td style="border: 1px solid #ccc; padding: 5px;"> twitter-2010 </td><td style="border: 1px solid #ccc; padding: 5px;"> 1.5B edges </td><td style="border: 1px solid #ccc; padding: 5px;"> Spark, 50 machines, 8.1 min </td><td style="border: 1px solid #ccc; padding: 5px;"> 13 min </td><td style="border: 1px solid #ccc; padding: 5px;"> 1 </td></tr> <tr><td style="border: 1px solid #ccc; padding: 5px;"> Pagerank - 100 iterations </td><td style="border: 1px solid #ccc; padding: 5px;"> uk-union </td><td style="border: 1px solid #ccc; padding: 5px;"> 3.8B edges </td><td style="border: 1px solid #ccc; padding: 5px;"> Stanford GPS (Pregel), 30 machines, 144 min </td><td style="border: 1px solid #ccc; padding: 5px;"> 581 min </td><td style="border: 1px solid #ccc; padding: 5px;"> 2  </td></tr> <tr><td style="border: 1px solid #ccc; padding: 5px;"> Web-graph Belief Propagation (1 iter.) </td><td style="border: 1px solid #ccc; padding: 5px;"> yahoo-web </td><td style="border: 1px solid #ccc; padding: 5px;"> 6.7B edges </td><td style="border: 1px solid #ccc; padding: 5px;"> Pegasus, 100 machines, 22 min </td><td style="border: 1px solid #ccc; padding: 5px;"> 27 min </td><td style="border: 1px solid #ccc; padding: 5px;"> 3 </td></tr> <tr><td style="border: 1px solid #ccc; padding: 5px;"> Matrix factorization (ALS), 10 iters </td><td style="border: 1px solid #ccc; padding: 5px;"> Netflix </td><td style="border: 1px solid #ccc; padding: 5px;"> 99M edges </td><td style="border: 1px solid #ccc; padding: 5px;"> GraphLab, 8-core machine, 4.7 min </td><td style="border: 1px solid #ccc; padding: 5px;"> 9.8 min </td><td style="border: 1px solid #ccc; padding: 5px;"> 4 </td></tr> <tr><td style="border: 1px solid #ccc; padding: 5px;"> Triangle counting </td><td style="border: 1px solid #ccc; padding: 5px;"> twitter-2010 </td><td style="border: 1px solid #ccc; padding: 5px;"> 1.5B edges </td><td style="border: 1px solid #ccc; padding: 5px;"> Hadoop, 1636 machines, 423 mins </td><td style="border: 1px solid #ccc; padding: 5px;"> 55 min </td><td style="border: 1px solid #ccc; padding: 5px;"> 5 </td></tr> </table>

**Performance on hard drive:** We repeated the same experiments on hard-drive, and the running times are roughly double compared to SSD. The hard-drive performance is thus sufficient for many purposes.

**Configuration:** 
     membudget_mb 3000  execthreads 2 loadthreads 2 niothreads 2 io.blocksize 1048576

<table>
  <tr><td>1</td><td>I. Stanton and G. Kliot. Streaming graph partitioning for large distributed graphs. 2012.</td></tr>
  <tr><td>2</td><td> S.Salihoglu and J.Widom. GPS: A Graph Processing System. Technical Report, pages 1–32, Apr. 2012.</td></tr>
  <tr><td>3</td><td> U. Kang, D. H. Chau, and C. Faloutsos. Inference of Beliefs on Billion-Scale Graphs. KDD-LDMTA’10, pages 1–7, June 2010. </td></tr>
  <tr><td>4</td><td> http://graphlab.org/datasets.html (Retrieved June 30, 2012) </td></tr>
  <tr><td>5</td><td> S. Suri and S. Vassilvitskii. Counting triangles and the curse of the last reducer. In Proceedings of the 20th international conference on World wide web, pages 607–614. ACM, 2011. </td></tr>
</table>

### Comparison to Giraph 

Apache [Giraph](http://giraph.apache.org/) is an open-source implementation of the Pregel graph engine, built on top of Hadoop. Based on a recent talk by the main developer of Giraph ( http://www.youtube.com/watch?v=b5Qmz4zPj-M ), Giraph running with 20 workers can run five iterations of !PageRank  on a graph with 5 billion edges in approx. 75 minutes. Estimated from our results above, GraphChi can execute similar task in roughly the same time on just one machine. The structure of the input graph affects the runtime, so a direct comparison is not possible without using the same input.  

### Preprocessing times

<table>
  <tr><td><b>Graph</b></td><td><b>Preprocessing time</b></td></tr>
  <tr><td>Netflix</td><td>1 min</td></tr>
  <tr><td>Twitter-2010</td><td>10 min</td></tr>
  <tr><td>uk-union</td><td>33 min</td></tr>
  <tr><td>yahoo-web</td><td>37 min</td></tr>
</table>
