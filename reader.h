// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef READER_H_
#define READER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <type_traits>
#include <sstream>
#include <vector>
#include <iterator>

#include "print_util.h"
#include "pvector.h"


/*
GAP Benchmark Suite
Class:  Reader
Author: Scott Beamer

Given filename, returns an edgelist or the entire graph (if serialized)
 - Intended to be called from Builder
 - Determines file format from the filename's suffix
 - If the input graph is serialized (.sg or .wsg), reads the graph
   directly into the returned graph instance
 - Otherwise, reads the file and returns an edgelist
*/


template <typename NodeID_, typename DestID_ = NodeID_,
          typename WeightT_ = NodeID_, bool invert = true>
class Reader {
  typedef EdgePair<NodeID_, DestID_> Edge;
  typedef pvector<Edge> EdgeList;
  std::string filename_;

 public:
  explicit Reader(std::string filename) : filename_(filename) {}

  std::string GetSuffix() {
    std::size_t suff_pos = filename_.rfind('.');
    if (suff_pos == std::string::npos) {
      std::cout << "Could't find suffix of " << filename_ << std::endl;
      std::exit(-1);
    }
    return filename_.substr(suff_pos);
  }


  EdgeList ReadInSNAP(std::ifstream &in) 
  {
    std::string line;
    std::vector<std::string> tokens;
    int num_vertices=0, num_edges=0;
    //proceed until you get out of '#' comments on the top
    do
    {
      getline(in, line);
      std::istringstream liness(line);
      tokens.clear();
      std::copy( std::istream_iterator<std::string>(liness),
          std::istream_iterator<std::string>(),
          back_inserter(tokens));
      if(tokens[1] == "Nodes:")
      {
        num_vertices = stoi(tokens[2]);
        num_edges    = stoi(tokens[4]);
      }
    }while(tokens[0].c_str()[0]=='#');

    //EdgeList el(num_edges);
    EdgeList el(num_edges);

    while (getline(in, line))
    {
      std::istringstream liness(line);
      tokens.clear();
      std::copy( std::istream_iterator<std::string>(liness),
          std::istream_iterator<std::string>(),
          back_inserter(tokens));

      uint64_t src = stoi(tokens[0]);
      uint64_t dst = stoi(tokens[0]);
      el.push_back(Edge(src, dst));
    }

    std::cout<<"edge size on the file : "<<num_edges<<" on the list : "<<el.size()<<std::endl;
    exit(0);
  }




  EdgeList ReadInEL(std::ifstream &in) {
    EdgeList el;
    NodeID_ u, v;
    while (in >> u >> v) {
      el.push_back(Edge(u, v));
    }
    return el;
  }

  EdgeList ReadInWEL(std::ifstream &in) {
    EdgeList el;
    NodeID_ u;
    NodeWeight<NodeID_, WeightT_> v;
    while (in >> u >> v) {
      el.push_back(Edge(u, v));
    }
    return el;
  }

  EdgeList ReadInGR(std::ifstream &in) {
    EdgeList el;
    char c;
    NodeID_ u;
    NodeWeight<NodeID_, WeightT_> v;
    while (!in.eof()) {
      c = in.peek();
      if (c == 'a') {
        in >> c >> u >> v;
        el.push_back(Edge(u, v));
      } else {
        in.ignore(200, '\n');
      }
    }
    return el;
  }

  EdgeList ReadFile(bool &needs_weights) {
    Timer t;
    t.Start();
    EdgeList el;
    std::string suffix = GetSuffix();
    std::ifstream file(filename_);
    if (!file.is_open()) {
      std::cout << "Couldn't open file " << filename_ << std::endl;
      std::exit(-2);
    }
    if (suffix == ".el" || suffix == ".net") {
      el = ReadInEL(file);
    } else if (suffix == ".wel") {
      needs_weights = false;
      el = ReadInWEL(file);
    } else if (suffix == ".gr") {
      needs_weights = false;
      el = ReadInGR(file);
    } else if (suffix == ".txt") {
      needs_weights = false;
      el = ReadInSNAP(file);
    } else {
      std::cout << "Unrecognized suffix: " << suffix << std::endl;
      std::exit(-3);
    }
    file.close();
    t.Stop();
    PrintTime("Read Time", t.Seconds());
    return el;
  }
/*
  CSRGraph<NodeID_, DestID_, invert> ReadSerializedGraph() {
    bool weighted = GetSuffix() == ".wsg";
    if (!std::is_same<NodeID_, SGID>::value) {
      std::cout << "serialized graphs only allowed for 32bit" << std::endl;
      std::exit(-5);
    }
    if (!weighted && !std::is_same<NodeID_, DestID_>::value) {
      std::cout << ".sg not allowed for weighted graphs" << std::endl;
      std::exit(-5);
    }
    if (weighted && std::is_same<NodeID_, DestID_>::value) {
      std::cout << ".wsg only allowed for weighted graphs" << std::endl;
      std::exit(-5);
    }
    if (weighted && !std::is_same<WeightT_, SGID>::value) {
      std::cout << ".wsg only allowed for int32_t weights" << std::endl;
      std::exit(-5);
    }
    std::ifstream file(filename_);
    if (!file.is_open()) {
      std::cout << "Couldn't open file " << filename_ << std::endl;
      std::exit(-6);
    }
    Timer t;
    t.Start();
    bool directed;
    SGOffset num_nodes, num_edges;
    DestID_ **index = nullptr, **inv_index = nullptr;
    DestID_ *neighs = nullptr, *inv_neighs = nullptr;
    file.read(reinterpret_cast<char*>(&directed), sizeof(bool));
    file.read(reinterpret_cast<char*>(&num_edges), sizeof(SGOffset));
    file.read(reinterpret_cast<char*>(&num_nodes), sizeof(SGOffset));
    pvector<SGOffset> offsets(num_nodes+1);
    neighs = new DestID_[num_edges];
    std::streamsize num_index_bytes = (num_nodes+1) * sizeof(SGOffset);
    std::streamsize num_neigh_bytes = num_edges * sizeof(DestID_);
    file.read(reinterpret_cast<char*>(offsets.data()), num_index_bytes);
    file.read(reinterpret_cast<char*>(neighs), num_neigh_bytes);
    index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
    if (directed && invert) {
      inv_neighs = new DestID_[num_edges];
      file.read(reinterpret_cast<char*>(offsets.data()), num_index_bytes);
      file.read(reinterpret_cast<char*>(inv_neighs), num_neigh_bytes);
      inv_index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, inv_neighs);
    }
    file.close();
    t.Stop();
    PrintTime("Read Time", t.Seconds());
    if (directed)
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes, index, neighs,
                                                inv_index, inv_neighs);
    else
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes, index, neighs);
  }
  */
};

#endif  // READER_H_
