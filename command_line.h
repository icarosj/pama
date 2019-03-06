// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef COMMAND_LINE_H_
#define COMMAND_LINE_H_

#include <getopt.h>

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <stdio.h>


/*
GAP Benchmark Suite
Class:  CLBase
Author: Scott Beamer

Handles command line argument parsing
 - Through inheritance, can add more options to object
 - For example, most kernels will use CLApp
*/


class CLBase {
 protected:
  int argc_;
  char** argv_;
  std::string name_;
  std::string get_args_ = "f:g:hsu:d:L";
  std::vector<std::string> help_strings_;

  int scale_ = -1;
  std::string filename_ = "";
  std::string dbpath_ = "./tmpdb";
  std::vector<std::string> files_;
  bool symmetrize_ = false;
  bool uniform_ = false;
  bool loadDB_ = false;

  void AddHelpLine(char opt, std::string opt_arg, std::string text,
                   std::string def = "") {
    const int kBufLen = 100;
    char buf[kBufLen];
    if (opt_arg != "")
      opt_arg = "<" + opt_arg + ">";
    if (def != "")
      def = "[" + def + "]";
    sprintf(buf, " -%c %-9s: %-57s%7s", opt, opt_arg.c_str(),
            text.c_str(), def.c_str());
    help_strings_.push_back(buf);
  }

 public:
  CLBase(int argc, char** argv, std::string name = "") :
         argc_(argc), argv_(argv), name_(name) {
    AddHelpLine('h', "", "print this help message");
    AddHelpLine('f', "file", "load graph from file. not available any more. use L)");
    AddHelpLine('s', "", "symmetrize input edge list", "false");
    AddHelpLine('g', "scale", "generate 2^scale kronecker graph");
    AddHelpLine('u', "scale", "generate 2^scale uniform-random graph");
    AddHelpLine('L', "load", "load files to form a snapshot");
    AddHelpLine('d', "db", "provide dir for mmap");
  }

  bool ParseArgs() {
    signed char c_opt;
    extern char *optarg;          // from and for getopt
    while ((c_opt = getopt(argc_, argv_, get_args_.c_str())) != -1) {
      HandleArg(c_opt, optarg);
    }

    for(int i=optind;i<argc_;i++)
    {
      files_.push_back(std::string(argv_[i]));
    }

    //if ((filename_ == "") && (scale_ == -1)) {
    if ((loadDB_==0) && (files_.size()== 0) && (filename_ == "") && (scale_ == -1)) {
      std::cout << "No graph input specified. (Use -h for help)" << std::endl;
      return false;
    }

    //make db directory 

    struct stat st;
    if(dbpath_ != "")
    {
      if (stat(dbpath_.c_str(), &st)) {
        if (mkdir(dbpath_.c_str(), 0755)) {
          abort();
        }
      }
    }


    if (scale_ != -1)
      symmetrize_ = true;
    return true;
  }

  void virtual HandleArg(signed char opt, char* opt_arg) {
    switch (opt) {
      case 'f': filename_ = std::string(opt_arg);           break;
      case 'd': dbpath_ = std::string(opt_arg);           break;
      case 'g': scale_ = atoi(opt_arg);                     break;
      case 'h': PrintUsage();                               break;
      case 's': symmetrize_ = true;                         break;
      case 'u': uniform_ = true; scale_ = atoi(opt_arg);    break;
      case 'L': loadDB_  = true;                              break;
    }
  }

  void PrintUsage() {
    std::cout << name_ << std::endl;
    // std::sort(help_strings_.begin(), help_strings_.end());
    for (std::string h : help_strings_)
      std::cout << h << std::endl;
    std::exit(0);
  }

  int scale() const { return scale_; }
  std::string filename() const { return filename_; }
  std::vector<std::string> files() const { return files_; }
  bool symmetrize() const { return symmetrize_; }
  bool uniform() const { return uniform_; }
  std::string dbpath() const {return dbpath_;}
  bool loadDB() const {return loadDB_;}
};



class CLApp : public CLBase {
  bool do_analysis_ = false;
  int num_trials_ = 1;
  int64_t start_vertex_ = -1;

 public:
  CLApp(int argc, char** argv, std::string name) : CLBase(argc, argv, name) {
    get_args_ += "an:r:";
    char buf[30];
    sprintf(buf, "%d", num_trials_);
    AddHelpLine('a', "a", "output analysis of last run", "false");
    AddHelpLine('n', "n", "perform n trials", buf);
    AddHelpLine('r', "node", "start from node r", "rand");
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'a': do_analysis_ = true;                    break;
      case 'n': num_trials_ = atoi(opt_arg);            break;
      case 'r': start_vertex_ = atol(opt_arg);          break;
      default: CLBase::HandleArg(opt, opt_arg);
    }
  }

  bool do_analysis() const { return do_analysis_; }
  int num_trials() const { return num_trials_; }
  int64_t start_vertex() const { return start_vertex_; }
};



class CLIterApp : public CLApp {
  int num_iters_;

 public:
  CLIterApp(int argc, char** argv, std::string name, int num_iters) :
    CLApp(argc, argv, name), num_iters_(num_iters) {
    get_args_ += "k:";
    char buf[15];
    sprintf(buf, "%d", num_iters_);
    AddHelpLine('k', "k", "perform k iterations", buf);
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'k': num_iters_ = atoi(opt_arg);            break;
      default: CLApp::HandleArg(opt, opt_arg);
    }
  }

  int num_iters() const { return num_iters_; }
};



class CLDelta : public CLApp {
  int delta_ = 1;

 public:
  CLDelta(int argc, char** argv, std::string name) : CLApp(argc, argv, name) {
    get_args_ += "d:";
    char buf[15];
    sprintf(buf, "%d", delta_);
    AddHelpLine('d', "d", "delta parameter", buf);
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'd': delta_ = atoi(opt_arg);               break;
      default: CLApp::HandleArg(opt, opt_arg);
    }
  }

  int delta() const { return delta_; }
};



class CLConvert : public CLBase {
  std::string out_filename_ = "";
  bool out_weighted_ = false;
  bool out_el_ = false;
  bool out_sg_ = false;

 public:
  CLConvert(int argc, char** argv, std::string name)
      : CLBase(argc, argv, name) {
    get_args_ += "e:b:w";
    AddHelpLine('b', "file", "output serialized graph to file");
    AddHelpLine('e', "file", "output edge list to file");
    AddHelpLine('w', "file", "make output weighted");
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'b': out_sg_ = true; out_filename_ = std::string(opt_arg);   break;
      case 'e': out_el_ = true; out_filename_ = std::string(opt_arg);   break;
      case 'w': out_weighted_ = true;                               break;
      default: CLBase::HandleArg(opt, opt_arg);
    }
  }

  std::string out_filename() const { return out_filename_; }
  bool out_weighted() const { return out_weighted_; }
  bool out_el() const { return out_el_; }
  bool out_sg() const { return out_sg_; }
};

#endif  // COMMAND_LINE_H_
