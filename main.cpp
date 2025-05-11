#include <igraph/igraph.h>
#include <algorithm>
#include <iostream>
#include <time.h>
#include <queue>
#include <vector>
#include <utility>
#include <ctime>
#include <cstdlib>
#include <random>
#include <stdio.h>
#include <libgen.h>
#include <unistd.h>
#include <fstream>

using namespace std;

int maxRepeating(igraph_vector_t& arr, int n, int k)
{
    // Iterate though input array, for every element
    // arr[i], increment arr[arr[i]%k] by k
    for (int i = 0; i< n; i++)
        VECTOR(arr)[int(VECTOR(arr)[i])%k] += k;

    // Find index of the maximum repeating element
    int max = VECTOR(arr)[0], result = 0;
    for (int i = 1; i < n; i++)
    {
        if (VECTOR(arr)[i] > max)
        {
            max = VECTOR(arr)[i];
            result = i;
        }
    }


    return result;
}





int simulate_thresholds(igraph_t &largest_c_graph, int max_comp_size, float gamma, int N, int E, int real_E, vector<int> & seed_vector)
{
  igraph_bool_t simple;

  igraph_is_simple(&largest_c_graph, &simple);
  if(simple == false)
  {
    cout<<"not simple graph, error\n";
    return -1;
  }


  std::priority_queue<std::pair<double, int> > q;
  igraph_vector_int_t result;
  igraph_vector_int_init(&result, 0);

  std::vector<int> initial_act(max_comp_size, 0);
  std::vector<int> number_of_decisions(max_comp_size, 0);
  std::vector<int> number_of_decisions_next(max_comp_size, 0);
  std::vector<int> initial_act_copy(max_comp_size, 0);
  std::vector<int> initial_act_copy2(max_comp_size, 0);

  //std::vector<std::pair<float, float> > thresholds = { {0.4,1.1}, {0.4,1.0}, {0.4,0.5}, {0.4,0.6} };
  //  std::vector<std::pair<float, float> > thresholds = { {0.5,1.1}, {0.5,1.0}, {0.5,0.8}, {0.2,0.6}, {0.5, 0.9} };
   std::vector<std::pair<float, float> > thresholds = { {0.5,1.1}, {0.2,1.1}};

  igraph_arpack_options_t arpack_options;
  igraph_arpack_options_init(&arpack_options);
  //igraph_pagerank(&largest_c_graph, IGRAPH_PAGERANK_ALGO_ARPACK, &result, 0, igraph_vss_all(),0,0.85,0,&arpack_options);


  igraph_degree(&largest_c_graph, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);

  for (int i = 0; i < igraph_vector_int_size(&result); ++i) {
    q.push(std::pair<double, int>(VECTOR(result)[i], i));
  }


  float thr1;
  float thr2;

  igraph_vector_int_t neighbors;
  igraph_vector_int_t neighbors2;

  igraph_vector_int_init(&neighbors, 0);
  igraph_vector_int_init(&neighbors2, 0);


  //graph Attributes
  igraph_real_t density;
  igraph_density(&largest_c_graph, &density, false);

  igraph_real_t transitivity;
  igraph_transitivity_undirected(&largest_c_graph, &transitivity, IGRAPH_TRANSITIVITY_ZERO);

  igraph_real_t assortativity_degree;
  igraph_assortativity_degree(&largest_c_graph, &assortativity_degree,1);

  std::vector<int> shuffle_nodes(max_comp_size, 0);
  for(int i = 0; i < max_comp_size; i++)
  shuffle_nodes[i]=i;
  mt19937 g(static_cast<uint32_t>(time(0)));
  //srand ( unsigned ( std::time(0) ) );
  shuffle ( shuffle_nodes.begin(), shuffle_nodes.end(),g );

  /*
  for(int i = 0; i < 10; i++)
  cout<<shuffle_nodes[i]<<" ";
  cout<<"\n";
  for(int i = 0; i < 10; i++)
  cout<<shuffle_nodes[max_comp_size-i-1]<<" ";
  cout<<"\n";
  return 0;
  */
for (int i = 0; i < seed_vector.size(); i++) {
    int ki = seed_vector[i];
    //std::cout << "index[" << ki << "] = " << "---"  << std::endl;

    initial_act[ki] = 1;
    initial_act_copy[ki] = 1;
    initial_act_copy2[ki] = 1;

}
  /*

  int k = 10; // number of indices we need, top 10 by degree
  for (int i = 0; i < k; ++i) {
    int ki = q.top().second;
    //ki = shuffle_nodes[i];
    //std::cout << "index[" << ki << "] = " << q.top().first  << std::endl;
    std::cout << "index[" << ki << "] = " << VECTOR(result)[ki]  << std::endl;

    initial_act[ki] = 1;
    initial_act_copy[ki] = 1;
    initial_act_copy2[ki] = 1;


    q.pop();
  }
*/

  //cout<<"path started\n";
  igraph_real_t avg_shortest_path;
  //igraph_average_path_length(&largest_c_graph, &avg_shortest_path,IGRAPH_UNDIRECTED, 1);
  //cout<<"path ready\n";

  igraph_matrix_t path_matrix;
  igraph_matrix_init(&path_matrix,0,0);

  igraph_matrix_t laplacian;
  igraph_matrix_init(&laplacian,0,0);

  igraph_vector_t from;
  igraph_vector_init(&from, 1000);
  igraph_vector_t to;
  igraph_vector_init(&to, 1000);
  for(int i = 0; i < 1000; i++)
  VECTOR(from)[i] = shuffle_nodes[i];
  for(int i = 0; i < 1000; i++)
  VECTOR(to)[i] = shuffle_nodes[max_comp_size-i-1];



  //igraph_shortest_paths(&largest_c_graph,&path_matrix,igraph_vss_vector(&from),igraph_vss_vector(&to),IGRAPH_ALL);

  float sample_avg = 0.0;
  int diam = 0;
/*
 for(int i = 0; i < 1000; i++)
  for(int j = 0; j < 1000; j++)
  {
  sample_avg+=MATRIX(path_matrix,i,j);
  if (MATRIX(path_matrix,i,j) > diam)
  diam = MATRIX(path_matrix,i,j);
  }
  sample_avg/=1000000;
  //cout<<avg_shortest_path<<" "<<sample_avg<<" "<<diam<<endl;
*/

  //igraph_laplacian(&largest_c_graph,&laplacian,NULL,false,NULL);

  igraph_vector_t eigenvalues;
  igraph_vector_init(&eigenvalues, 0);

  //igraph_lapack_dsyevr(&laplacian, IGRAPH_LAPACK_DSYEV_ALL, 0,0,0,0,0,1e-10,&eigenvalues,NULL,NULL);


  //for(int i = 0; i < igraph_vector_size(&eigenvalues); i++)
  {
    //std::cout<<i<<" "<<VECTOR(eigenvalues)[i]<<std::endl;
  }

  //std::cout<<"start\n";
  //std::cout<<VECTOR(eigenvalues)[1]<<std::endl;
  //return 0;
  //one decision per node, initial idea
  for(int run = 0; run < thresholds.size(); run++)
  {
   thr1 = thresholds[run].first;
   thr2 = thresholds[run].second;

   initial_act = initial_act_copy2;
   initial_act_copy = initial_act_copy2;

   std::vector<float> res_line ={};
   res_line.push_back(1);
   res_line.push_back(N);
   res_line.push_back(E);
   res_line.push_back(max_comp_size);
   res_line.push_back(real_E);

//   res_line.push_back(VECTOR(eigenvalues)[1]);
  res_line.push_back(0);

   res_line.push_back(sample_avg);
   res_line.push_back(diam);

   res_line.push_back(density);
   res_line.push_back(transitivity);
   res_line.push_back(assortativity_degree);

   res_line.push_back(gamma);
   res_line.push_back(thr1);
   res_line.push_back(thr2);
  //  res_line.push_back(i);
//  cout<<"start\n";
  for(int i = 0; i < 250; i++)
  {
    std::vector<int> decided(max_comp_size, 0);
//cout<<"iter\n";

      for(int seed = 0; seed < max_comp_size; seed++)
      {
        if (initial_act[seed] != 1)
        continue;

        igraph_neighbors(&largest_c_graph, &neighbors, seed, IGRAPH_ALL);
        for(int j = 0; j < igraph_vector_int_size(&neighbors); j++)
        {
          if (decided[ VECTOR(neighbors)[j]] != 0)
          continue;
          decided[ VECTOR(neighbors)[j]] =1;
          igraph_neighbors(&largest_c_graph, &neighbors2, VECTOR(neighbors)[j], IGRAPH_ALL);
          int cnt_act = 0;
          int cnt_nonact = 0;
          for(int k = 0; k < igraph_vector_int_size(&neighbors2); k++)
          {

            if (initial_act[VECTOR(neighbors2)[k]] == 1)
            cnt_act++;
            else
            cnt_nonact++;
          }

          float pt = float(cnt_act)/(cnt_act+cnt_nonact);
          //std::cout<<pt<<std::endl;
          if (initial_act[VECTOR(neighbors)[j]] == 0)
          {
            //std::cout<<"here";
            if (pt >= thr1 && pt < thr2)
            initial_act_copy[VECTOR(neighbors)[j]] = 1;
            else if (pt >= thr2)
            initial_act_copy[VECTOR(neighbors)[j]] = -1;
          }

        }


      }
    if (initial_act == initial_act_copy)
    {

      //std::cout<<"equal "<<i<<std::endl;
      int cnt = 0;
      for(int i = 0;i < initial_act.size(); i ++)
      if (initial_act[i] == 1)
      cnt++;
      //std::cout<<"step "<<i<<" thr1 = "<<thr1<<" thr2 = "<<thr2<<" activated = "<<cnt<<std::endl;

      res_line.push_back(cnt);

      break;
    }
    int cnt = 0;
    for(int i = 0;i < initial_act.size(); i ++)
    if (initial_act[i] == 1)
    cnt++;
    //std::cout<<"step "<<i<<" thr1 = "<<thr1<<" thr2 = "<<thr2<<" activated = "<<cnt<<std::endl;
    res_line.push_back(cnt);

    initial_act = initial_act_copy;


  }
  for (auto i: res_line)
   std::cout << i << ' ';
   std::cout << std::endl;
 }


 //more than one decision
 for(int run = 0; run < thresholds.size(); run++)
 {
  thr1 = thresholds[run].first;
  thr2 = thresholds[run].second;

  initial_act = initial_act_copy2;
  initial_act_copy = initial_act_copy2;
  fill(number_of_decisions.begin(),number_of_decisions.end(),0);
  fill(number_of_decisions_next.begin(),number_of_decisions_next.end(),0);

  std::vector<float> res_line ={};
  res_line.push_back(2);
  res_line.push_back(N);
  res_line.push_back(E);
  res_line.push_back(max_comp_size);
  res_line.push_back(real_E);

  //res_line.push_back(VECTOR(eigenvalues)[1]);
  res_line.push_back(0);

  res_line.push_back(sample_avg);
  res_line.push_back(diam);

  res_line.push_back(density);
  res_line.push_back(transitivity);
  res_line.push_back(assortativity_degree);


  res_line.push_back(gamma);
  res_line.push_back(thr1);
  res_line.push_back(thr2);
 //  res_line.push_back(i);

 for(int i = 0; i < 250; i++)
 {
std::vector<int> decided(max_comp_size, 0);

     for(int seed = 0; seed < max_comp_size; seed++)
     {
       if (initial_act[seed] != 1)
       continue;

       igraph_neighbors(&largest_c_graph, &neighbors, seed, IGRAPH_ALL);
       for(int j = 0; j < igraph_vector_int_size(&neighbors); j++)
       {
         if (decided[ VECTOR(neighbors)[j]] != 0)
         continue;
         decided[ VECTOR(neighbors)[j]] =1;

         igraph_neighbors(&largest_c_graph, &neighbors2, VECTOR(neighbors)[j], IGRAPH_ALL);
         int cnt_act = 0;
         int cnt_nonact = 0;
         for(int k = 0; k < igraph_vector_int_size(&neighbors2); k++)
         {

           if (initial_act[VECTOR(neighbors2)[k]] == 1)
           cnt_act++;
           else
           cnt_nonact++;
         }

         float pt = float(cnt_act)/(cnt_act+cnt_nonact);
         //std::cout<<pt<<std::endl;

         if (number_of_decisions[VECTOR(neighbors)[j]] <= 1)
         {

           //if (initial_act[VECTOR(neighbors)[j]] == 0)
           {
             if (pt >= thr1 && pt < thr2)
             initial_act_copy[VECTOR(neighbors)[j]] = 1;
             else if (pt >= thr2)
             initial_act_copy[VECTOR(neighbors)[j]] = 0;

             if (initial_act[VECTOR(neighbors)[j]] != initial_act_copy[VECTOR(neighbors)[j]])
             number_of_decisions_next[VECTOR(neighbors)[j]] = 1+number_of_decisions[VECTOR(neighbors)[j]];
           }
           //else
           {


           }
         }



       }


     }
   if (initial_act == initial_act_copy)
   {

     //std::cout<<"equal "<<i<<std::endl;
     int cnt = 0;
     for(int i = 0;i < initial_act.size(); i ++)
     if (initial_act[i] == 1)
     cnt++;
     //std::cout<<"step "<<i<<" thr1 = "<<thr1<<" thr2 = "<<thr2<<" activated = "<<cnt<<std::endl;

     res_line.push_back(cnt);

     break;
   }
   int cnt = 0;
   for(int i = 0;i < initial_act.size(); i ++)
   if (initial_act[i] == 1)
   cnt++;
   //std::cout<<"step "<<i<<" thr1 = "<<thr1<<" thr2 = "<<thr2<<" activated = "<<cnt<<std::endl;
   res_line.push_back(cnt);

   initial_act = initial_act_copy;
   number_of_decisions = number_of_decisions_next;

 }
 for (auto i: res_line)
  std::cout << i << ' ';
  std::cout << std::endl;
}


  return 0;
}

int main(int argc, char *argv[])
{

      float gamma = stof(argv[1]);
      char* graph_file_name = argv[2];
      char* seed_file_name = argv[3];

      vector<int> seed_vector;

      // read seed file
      string line;
  ifstream myfile (seed_file_name);
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      //cout << line << '\n';
      seed_vector.push_back(stoi(line));
    }
    myfile.close();
  }

  else cout << "Unable to open file"; 
     
      // end seed file reading


     std::clock_t begin = std::clock();

  //  for (float gamma = 2.0;  gamma <= 3.6; gamma+=0.05)
{
     igraph_integer_t diameter;
     igraph_integer_t num_comp;
     igraph_t graph;
     igraph_t largest_c_graph;
     igraph_vector_int_t vertices;
     igraph_vector_int_t components;
     igraph_vector_int_t new_vertices;

     igraph_vector_int_init(&vertices, 0);
     igraph_vector_int_init(&components, 0);

     int N = 1000000;
     int E = 100*N;
     //float gamma = 2;
     int exponent_in = -1;

     igraph_rng_seed(igraph_rng_default(),  time(0));


     //igraph_static_power_law_game(&graph, N, E, gamma, exponent_in, false,false,false);


     FILE *ifile;
     ifile = fopen(graph_file_name,"r");
     //igraph_read_graph_pajek(&graph,ifile);
     igraph_read_graph_edgelist(&graph,ifile,(int)gamma,false);// 15233
     //std::cout<<"read"<<std::endl;

     igraph_bool_t simple;
     //cout<<"graph done\n";
     igraph_is_simple(&graph, &simple);
     if(simple == false)
     {
       cout<<"not simple graph, simplifying\n";
       igraph_simplify(&graph,1,1,0);
     }

     igraph_clusters(&graph, &vertices, &components, &num_comp,IGRAPH_WEAK);
     //std::cout<<num_comp<<std::endl;

     /*
     for(int i = 0; i < igraph_vector_size(&components); i++)
     {
       std::cout<<i<<" "<<VECTOR(components)[i]<<std::endl;
     }
     */

    int max_comp_size = (int)igraph_vector_int_max(&components);
    int max_comp_num = (int)igraph_vector_int_which_max(&components);
    //std::cout<<max_comp_size<<" "<<max_comp_num<<std::endl;
    igraph_vector_int_init(&new_vertices, max_comp_size);



     int j = 0;
     for(int i = 0; i < igraph_vector_int_size(&vertices); i++)
     {

       //std::cout<<i<<" "<<VECTOR(vertices)[i]<<std::endl;
       if (VECTOR(vertices)[i] == max_comp_num)
          VECTOR(new_vertices)[j++] = i;
     }

     igraph_induced_subgraph(&graph, &largest_c_graph, igraph_vss_vector(&new_vertices), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH);
     igraph_destroy(&graph);

     //std::cout<<igraph_vcount(&largest_c_graph)<<" "<<igraph_ecount(&largest_c_graph)<<std::endl;

     //cout<<"start sim\n";
     simulate_thresholds(largest_c_graph, max_comp_size, gamma, N, E, igraph_ecount(&largest_c_graph), seed_vector);

     //seeds






     igraph_vector_int_destroy(&vertices);
     igraph_vector_int_destroy(&components);
     igraph_destroy(&largest_c_graph);
   }
     std::clock_t end = std::clock();
     double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
     std::cout<<"total: "<<elapsed_secs<<" seconds\n";
     return 0;
}

