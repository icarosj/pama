//Re-implementing llama based on gapb


#ifndef SNAPSHOTS_H_
#define SNAPSHOTS_H_

#include "deltagraph.h"
#include "pvector.h"
#include <map>
#include <set>
#include <list>


//CONT record size divided by edgelist entry
#define CONT_SIZE (sizeof(delta_vt_elem)/sizeof(DestID_))


//a global function? are you sure??


//template <class NodeID_, class DestID_ = NodeID_, bool MakeInverse = true>
template <class NodeID_, class DestID_ , bool MakeInverse >
class Snapshots
{
  //copied from builder.h
  typedef EdgePair<NodeID_, DestID_> Edge;
  typedef pvector<Edge> EdgeList;
  private:
    int64_t newest_level=-1;


    //indirection : level_n_index to vertices.
    //i could also have direct pointers to vertices,
    //but i would need same structure for base and delta vertices
    //
    //don't need indirection for the base graph
    //
    //indirections[level][vid >> perpagebits]
    //HACK::always access it with level-1 (no indirection at the base graph)
    std::vector<std::vector<delta_vt_elem*>> out_indirections_;
    std::vector<std::vector<delta_vt_elem*>> in_indirections_;

    //back to deltas only
    std::vector<DeltaGraph<NodeID_>> deltas_;
    //DeltaGraph<NodeID_>* deltas_;
  public:
    //TODO: i need constructors

    uint64_t num_nodes(int level=-1) const
    {
      if(level==-1) level = newest_level;


      return deltas_[level].num_nodes();
    }


    uint64_t num_edges(int level=-1) const
    {
      if(level==-1) level = newest_level;
      return deltas_[level].num_edges();
    }

    void write_meta_file(std::string dbpath)
    {
      std::string filename = dbpath + "/meta";
      //write from the start
      std::ofstream of(filename, std::ofstream::out );

      int num_levels = newest_level+1;
      of<<num_levels<<std::endl<<std::endl;

      for(int i=0;i<num_levels;i++)
      {
        of<<deltas_[i].directed()<<std::endl;
        of<<deltas_[i].num_nodes()<<std::endl;
        of<<deltas_[i].num_edges()<<std::endl;
        of<<deltas_[i].out_delta_size()<<std::endl;
        of<<deltas_[i].in_delta_size()<<std::endl;
        of<<deltas_[i].out_et_size()<<std::endl;
        of<<deltas_[i].in_et_size()<<std::endl<<std::endl;
      }
      of.close();
    }


    void write_indirection_file(std::string dbpath)
    {

      int num_levels = newest_level+1;

      //level that the indirection is stored
      for(int v_level=0;v_level<num_levels;v_level++)
      {
        debug_print("writing indirection level %d\n", v_level);
        std::string out_filename = dbpath + "/out_indirec"
          + std::to_string(v_level);
        std::string in_filename = dbpath + "/in_indirec"
          + std::to_string(v_level);
        std::ofstream out_f(out_filename, std::ofstream::out );
        std::ofstream  in_f(in_filename, std::ofstream::out );

        uint64_t num_nodes = this->num_nodes(v_level);
        uint64_t num_pages = 1 + IND_V_PG(num_nodes-1);
        debug_print("num nodes %lu num pages %lu\n", num_nodes, num_pages);

        for(uint64_t v_pg=0;v_pg<num_pages;v_pg++)
        {
          delta_vt_elem* out_ind_entry = out_indirections_[v_level][v_pg];
          delta_vt_elem*  in_ind_entry =  in_indirections_[v_level][v_pg];

          //look through physical v tables to find level and index
          bool out_found=false, in_found=false;

          for(int p_level=0;p_level<num_levels;p_level++)
          {
            delta_vt_elem* out_delta =
              deltas_[p_level].out_delta_vt(0);
            delta_vt_elem* in_delta =
              deltas_[p_level].in_delta_vt(0);

            uint64_t  out_delta_size =
              deltas_[p_level].out_delta_size();
            uint64_t   in_delta_size =
              deltas_[p_level].in_delta_size();

            if( out_delta <= out_ind_entry &&
                out_ind_entry < out_delta + out_delta_size)
            { //found it
              level_n_index location = level_n_index(
                  p_level, out_ind_entry - out_delta);
              int64_t content = location.content();
              out_f.write((char*)&content, sizeof(content));
              out_found=1;
            }

            if( in_delta <= in_ind_entry &&
                in_ind_entry < in_delta + in_delta_size)
            { //found it
              level_n_index location = level_n_index(
                  p_level, in_ind_entry - in_delta);
              int64_t content = location.content();
              in_f.write((char*)&content, sizeof(content));
              in_found=1;
            }
          }
          ASSERT(out_found, ==, true);
          ASSERT(in_found, ==, true);

        }//for v_pg
        out_f.close();
        in_f.close();
       /* std::ifstream tmp(in_filename, std::ifstream::in);
        std::cout<<"testing file "<< out_filename<<std::endl;
        while(!tmp.eof())
        {
          int64_t content;
          tmp.read((char*)&content, sizeof(content));

          level_n_index tmplni = level_n_index(content);
          printf("level %u index %lu\n", tmplni.level(), tmplni.index());
        }
        tmp.close();
*/

      }//for v_level


    }//write_indirection_file

    int read_database(std::string dbpath)
    {
      std::string filename = dbpath + "/meta";
      //write from the start
      std::ifstream inmeta(filename, std::ofstream::in );

      int num_levels ;
      inmeta>>num_levels;

      for(int level=0;level<num_levels;level++)
      {
        bool directed;
        uint64_t num_nodes, num_edges,
                 out_delta_size, in_delta_size,
                 out_et_size, in_et_size;

        inmeta>>directed;
        inmeta>>num_nodes;
        inmeta>>num_edges;
        inmeta>>out_delta_size;
        inmeta>>in_delta_size;
        inmeta>>out_et_size;
        inmeta>>in_et_size;

        gen_delta_from_db(dbpath, level, directed,
            num_nodes, num_edges,
            out_delta_size, in_delta_size,
            out_et_size, in_et_size);

        read_indirection_file(dbpath, level, num_nodes, 1);
        read_indirection_file(dbpath, level, num_nodes, 0);
      }
      inmeta.close();
      newest_level = num_levels-1;
      return num_levels;
    }



    //TODO move to private
    void gen_delta_from_db(std::string dbpath,
        int level, bool directed,
        uint64_t num_nodes, uint64_t num_edges,
        uint64_t out_delta_size, uint64_t in_delta_size,
        uint64_t out_et_size, uint64_t in_et_size )
    {
      delta_vt_elem* out_delta_vt =
        (delta_vt_elem*)mmap_open("out_V_", dbpath, level,
          sizeof(delta_vt_elem) * out_delta_size);

      delta_vt_elem* in_delta_vt =
        (delta_vt_elem*)mmap_open("in_V_", dbpath, level,
            sizeof(delta_vt_elem) * in_delta_size);


      DestID_* out_et =
        (DestID_*)mmap_open("out_E_", dbpath, level,
            sizeof(DestID_) * out_et_size);
      DestID_* in_et =
        (DestID_*)mmap_open("in_E_", dbpath, level,
            sizeof(DestID_) *in_et_size);

      ASSERT(deltas_.size(), ==, level);

      DeltaGraph<NodeID_> dg
        (num_nodes,
         num_edges,
         out_delta_vt, in_delta_vt,
         out_et,  in_et,
         out_delta_size, in_delta_size,
         out_et_size, in_et_size
        );
      deltas_.push_back(dg);


    }


    void read_indirection_file(std::string dbpath, int level, int num_nodes, bool is_out)
    {
      std::string inout = is_out? "/out_indirec":"/in_indirec";
      std::string filename = dbpath + inout + std::to_string(level);
      std::ifstream ifs(filename, std::ifstream::in);

      uint64_t num_pages = 1 + IND_V_PG(num_nodes-1);
      std::cout<<"reading file "<< filename<<" num_pages "<<num_pages<<std::endl;
      std::vector<delta_vt_elem*> indirection;
      indirection.resize(num_pages,0);

      uint64_t v_pg=0;
      for(int i=0;i<num_pages;i++)
      {
        int64_t content;
        ifs.read((char*)&content, sizeof(content));

        level_n_index lni = level_n_index(content);
        uint32_t p_level = lni.level();
        uint64_t p_index = lni.index();

        if(is_out)
        { //out
          indirection[v_pg] = deltas_[p_level].out_delta_vt(p_index);
        }
        else
        {//in
          indirection[v_pg] = deltas_[p_level].in_delta_vt(p_index);
        }
        v_pg++;
      }


      if(is_out)
      {
        ASSERT(out_indirections_.size(), == , level);
        out_indirections_.push_back(indirection);
      }
      else
      {
        ASSERT(in_indirections_.size(), == , level);
        in_indirections_.push_back(indirection);
      }



      ifs.close();


    }//read_indirection_file











    void create_base_from_EL(EdgeList& el, std::string dbpath)
    {
      uint64_t next_num_nodes = 0;
      ASSERT(newest_level,==,-1);

      std::vector<delta_vt_elem*> out_indirection;
      std::vector<delta_vt_elem*> in_indirection;
/****STEP 0. follow the EL and count.**************/
      //can I parallelize it?

      std::set<uint64_t> new_out_pages;
      std::set<uint64_t> new_in_pages;
      std::map<NodeID_,std::list<DestID_>> out_edges;
      std::map<NodeID_,std::list<DestID_>> in_edges;

      int next_num_edges=count_edges(el, next_num_nodes, new_out_pages, new_in_pages,
          out_edges, in_edges);

/*****STEP 1.5. identify empty pages*************/
      //if vertices with large number ahead gets inserted,
      //we have holes the middle
      std::set<uint64_t> empty_out_pages;
      std::set<uint64_t> empty_in_pages;
      uint64_t next_pages= 1 + IND_V_PG(next_num_nodes-1);
      std::cout<<"we need "<<next_pages<<" pages"<<std::endl;

      ASSERT(next_num_nodes, !=, 0); //don't support empty deltas

      next_pages = count_empty_pages(0, next_pages,
          new_out_pages, new_in_pages,
          empty_out_pages, empty_in_pages);


      out_indirection.resize(next_pages,0);
      in_indirection.resize(next_pages,0);



/*****STEP 2. alloc the structures*******************/
      //input : new_io_pages  dbpath next_num_edges num_cont_record
      //output: io_delta_size(dbg)  io_delta_vt  io_et io_nullcont_idx

      delta_vt_elem* out_delta_vt, *in_delta_vt;
      DestID_* out_et , *in_et;

      uint64_t out_delta_size, in_delta_size;
      uint64_t out_et_size, in_et_size;

      uint64_t in_nullcont_idx ;
      uint64_t out_nullcont_idx;
      alloc_tables(
          new_out_pages, new_in_pages,
          empty_out_pages, empty_in_pages,
          out_edges, in_edges,
          next_num_edges, dbpath,
          true,//it is a base
          &out_delta_vt, &in_delta_vt, &out_et, &in_et,
          &out_delta_size, &in_delta_size,
          &out_et_size, &in_et_size,
          &in_nullcont_idx, &out_nullcont_idx);



      //set the indirections it for the base
      int i=0;
      for(auto iter=out_indirection.begin();
          iter!=out_indirection.end();i++,iter++)
      {
        (*iter) = out_delta_vt + i*IND_VERTICES_PER_PAGE;
      }
      i=0;
      for(auto iter=in_indirection.begin();
          iter!=in_indirection.end();i++,iter++)
      {
        (*iter) = in_delta_vt + i*IND_VERTICES_PER_PAGE;
      }



      /*****STEP 3.fill the empty VTs********/
      //virtual pages --indirec---> physical pages
      uint64_t p_pg = 0;
      {

        //init the pages. init each vt entry with null conts
        for(uint64_t v_pg: new_out_pages)
        {
          {
            uint64_t p_vid = p_pg*IND_VERTICES_PER_PAGE;
            out_indirection[v_pg] = out_delta_vt + p_pg*IND_VERTICES_PER_PAGE;
            for(int i=0;i<IND_VERTICES_PER_PAGE ;
                i++,
                p_vid = p_pg*IND_VERTICES_PER_PAGE + i)
            { //fill with null vts
              out_delta_vt[p_vid] =
                delta_vt_elem(newest_level+1, out_nullcont_idx, 0, 0);
            }
          }
          p_pg++;
        }

        //put empty pages
        if(empty_out_pages.size()!=0)
        {
          //create the empty page shared instance
          uint64_t p_vid = p_pg*IND_VERTICES_PER_PAGE;
          for(int i=0;i<IND_VERTICES_PER_PAGE ;
              i++,
              p_vid = p_pg*IND_VERTICES_PER_PAGE + i)
          { //fill with null vts
            out_delta_vt[p_vid] =
              delta_vt_elem(newest_level+1, out_nullcont_idx, 0, 0);
          }

          for(uint64_t v_pg: empty_out_pages)
          { //all to the shared empty page
            out_indirection[v_pg] = out_delta_vt +  p_pg*IND_VERTICES_PER_PAGE;
          }
          p_pg++;
        }
        ASSERT((p_pg*IND_VERTICES_PER_PAGE) ,==, out_delta_size);
      }



      //now do the same for in_vt
      p_pg=0;
      {

        //init the pages. init each vt entry with null conts
        for(uint64_t v_pg: new_in_pages)
        {
          {
            uint64_t p_vid = p_pg*IND_VERTICES_PER_PAGE;
            in_indirection[v_pg] = in_delta_vt + p_pg*IND_VERTICES_PER_PAGE;
            for(int i=0;i<IND_VERTICES_PER_PAGE ;
                i++,
                p_vid = p_pg*IND_VERTICES_PER_PAGE + i)
            { //fill with null vts
              in_delta_vt[p_vid] =
                delta_vt_elem(newest_level+1, in_nullcont_idx, 0, 0);
            }
          }
          p_pg++;
        }

        //put empty pages
        if(empty_in_pages.size()!=0)
        {
          //create the empty page shared instance
          uint64_t p_vid = p_pg*IND_VERTICES_PER_PAGE;
          for(int i=0;i<IND_VERTICES_PER_PAGE ;
              i++,
              p_vid = p_pg*IND_VERTICES_PER_PAGE + i)
          { //fill with null vts
            in_delta_vt[p_vid] =
              delta_vt_elem(newest_level+1, in_nullcont_idx, 0, 0);
          }

          for(uint64_t v_pg: empty_in_pages)
          { //all to the shared empty page
            in_indirection[v_pg] = in_delta_vt +  p_pg*IND_VERTICES_PER_PAGE;
          }
          p_pg++;
        }
        ASSERT((p_pg*IND_VERTICES_PER_PAGE) ,==, in_delta_size);
      }


      /******STEP 4. put the new edges*****/

      {
        uint64_t edge_index=0;
        for(auto elist : out_edges)
        {
          uint64_t start_edge_index=edge_index;
          NodeID_ src = elist.first;
          for(DestID_ dest : elist.second)
          {
            out_et[edge_index++] = dest;
          }

          //now replace the vt entry to the new edgetable
          delta_vt_elem* src_vt =
            out_indirection[IND_V_PG(src)]
            + (IND_V_OFFSET(src));

          *src_vt =
            delta_vt_elem(newest_level+1, start_edge_index,
                elist.second.size(), elist.second.size());

          ASSERT(elist.second.size() ,==,
              edge_index - start_edge_index);
        }


        ASSERT(edge_index, ==,out_nullcont_idx);
        //TODO there must be a way to eliminate this
        //nullcont at the base graph..
        *(delta_vt_elem*)(out_et + out_nullcont_idx) =
          CONT_NULL;
      }



      //now the in edges
      {
        uint64_t edge_index=0;
        for(auto elist : in_edges)
        {
          uint64_t start_edge_index=edge_index;
          NodeID_ src = elist.first;
          for(DestID_ dest : elist.second)
          {
            in_et[edge_index++] = dest;
          }

          //now replace the vt entry to the new edgetable
          delta_vt_elem* src_vt =
            in_indirection[IND_V_PG(src)]
            + (IND_V_OFFSET(src));

          *src_vt =
            delta_vt_elem(newest_level+1, start_edge_index,
                elist.second.size(), elist.second.size());

          ASSERT(elist.second.size() ,==,
              edge_index - start_edge_index);
        }


        ASSERT(edge_index, ==,in_nullcont_idx);
        //TODO there must be a way to eliminate this
        //nullcont at the base graph..
        *(delta_vt_elem*)(in_et + in_nullcont_idx) =
          CONT_NULL;
      }

/****STEP 5. add the indirection to the list************/
      out_indirections_.push_back(out_indirection);
      in_indirections_.push_back(in_indirection);


      DeltaGraph<NodeID_> dg
        (next_num_nodes,
         el.size(),
         out_delta_vt, in_delta_vt,
         out_et, in_et,
         out_delta_size, in_delta_size,
         out_et_size, in_et_size
         );
      deltas_.push_back(dg);


      newest_level++;
    }


    //and this one should be called on adding deltas only
    //TODO: split them into multiple functions
    void add_delta_from_EL(EdgeList& el, std::string dbpath)
    {
      uint64_t prev_num_nodes = num_nodes(newest_level);
      uint64_t next_num_nodes = prev_num_nodes;
      ASSERT(newest_level,!=,-1);
/****STEP 0. copy the indirection arrays*********/
      std::vector<delta_vt_elem*> out_indirection;
      std::vector<delta_vt_elem*> in_indirection;

      out_indirection = std::vector<delta_vt_elem*>(out_indirections_[newest_level]);
      in_indirection =  std::vector<delta_vt_elem*>(in_indirections_[newest_level]);

/****STEP 1. follow the EL and count.**************/
      //can I parallelize it?

      std::set<uint64_t> new_out_pages;
      std::set<uint64_t> new_in_pages;
      std::map<NodeID_,std::list<DestID_>> out_edges;
      std::map<NodeID_,std::list<DestID_>> in_edges;

      int next_num_edges=count_edges(el, next_num_nodes, new_out_pages, new_in_pages,
          out_edges, in_edges);



/*****STEP 1.5. identify empty pages*************/
      //if vertices with large number ahead gets inserted,
      //we have holes the middle
      std::set<uint64_t> empty_out_pages;
      std::set<uint64_t> empty_in_pages;

      uint64_t prev_pages= 1 + IND_V_PG(prev_num_nodes-1);
      uint64_t next_pages= 1 + IND_V_PG(next_num_nodes-1);

      ASSERT( prev_pages ,==,
          out_indirection.size());

      next_pages = count_empty_pages(prev_pages, next_pages,
          new_out_pages, new_in_pages,
          empty_out_pages, empty_in_pages);

      out_indirection.resize(next_pages,0);
      in_indirection.resize(next_pages,0);




/*****STEP 2. alloc the structures*******************/
      //input : new_io_pages  dbpath next_num_edges is_base
      //output: io_delta_size(dbg)  io_delta_vt  io_et io_nullcont_idx

      delta_vt_elem *out_delta_vt, *in_delta_vt;
      DestID_ *out_et , *in_et;

      uint64_t out_delta_size, in_delta_size;
      uint64_t out_et_size, in_et_size;

      uint64_t in_nullcont_idx ;
      uint64_t out_nullcont_idx;
      alloc_tables(
          new_out_pages, new_in_pages,
          empty_out_pages, empty_in_pages,
          out_edges, in_edges,
          next_num_edges, dbpath,
          false, //it is a delta. we need cont.records
          &out_delta_vt, &in_delta_vt, &out_et, &in_et,
          &out_delta_size,&in_delta_size,
          &out_et_size, &in_et_size,
          &in_nullcont_idx, &out_nullcont_idx);





      /*****STEP 3.COW the necessary VTs********/
      //virtual pages --indirec---> physical pages
      {
        uint64_t p_pg = 0;
        for(uint64_t v_pg: new_out_pages)
        {
          if(v_pg < prev_pages)
          { //existing pages

            delta_vt_elem* orig_delta_vt = out_indirection[v_pg];


            uint64_t v_vid = v_pg*IND_VERTICES_PER_PAGE;
            uint64_t p_vid = p_pg*IND_VERTICES_PER_PAGE;
            for(int i=0;i<IND_VERTICES_PER_PAGE ;
                //              && p_vid<num_nodes() ;
                i++,
                v_vid = v_pg*IND_VERTICES_PER_PAGE + i,
                p_vid = p_pg*IND_VERTICES_PER_PAGE + i)
            {
              if(v_vid < prev_num_nodes)
              {//normal copy
                out_delta_vt[p_vid] = orig_delta_vt[i];
              }
              else
              {//fill the rest with null vts
                out_delta_vt[p_vid] =
                  delta_vt_elem(newest_level+1, out_nullcont_idx, 0, 0);
              }
            }
          } //v_pg<prev_pages
          else //new pages. indirections are set to zero
          {
            uint64_t v_vid = v_pg*IND_VERTICES_PER_PAGE;
            uint64_t p_vid = p_pg*IND_VERTICES_PER_PAGE;
            for(int i=0;i<IND_VERTICES_PER_PAGE ;
                i++,
                v_vid = v_pg*IND_VERTICES_PER_PAGE + i,
                p_vid = p_pg*IND_VERTICES_PER_PAGE + i)
            { //fill with null vts
              ASSERT(v_vid ,>=, prev_num_nodes);
              out_delta_vt[p_vid] =
                delta_vt_elem(newest_level+1, out_nullcont_idx, 0, 0);
            }
          }
          //vt cow done. not re-point
          out_indirection[v_pg] =  out_delta_vt + p_pg*IND_VERTICES_PER_PAGE;
          p_pg++;
        }

        //put empty pages
        if(empty_out_pages.size()!=0)
        {
          //create the empty page shared instance
          uint64_t p_vid = p_pg*IND_VERTICES_PER_PAGE;
          for(int i=0;i<IND_VERTICES_PER_PAGE ;
              i++,
              p_vid = p_pg*IND_VERTICES_PER_PAGE + i)
          { //fill with null vts
            out_delta_vt[p_vid] =
              delta_vt_elem(newest_level+1, out_nullcont_idx, 0, 0);
          }

          for(uint64_t v_pg: empty_out_pages)
          { //all to the shared empty page
            ASSERT(v_pg*IND_VERTICES_PER_PAGE ,>=, prev_num_nodes);
            out_indirection[v_pg] = out_delta_vt + p_pg*IND_VERTICES_PER_PAGE;
          }
          p_pg++;
        }
        ASSERT((p_pg*IND_VERTICES_PER_PAGE) ,==, out_delta_size);
      }



      //now do the same for in_vt
      {
        uint64_t p_pg = 0;
        for(uint64_t v_pg: new_in_pages)
        {
          if(v_pg < prev_pages)
          { //existing pages

            delta_vt_elem* orig_delta_vt = in_indirection[v_pg];


            uint64_t v_vid = v_pg*IND_VERTICES_PER_PAGE;
            uint64_t p_vid = p_pg*IND_VERTICES_PER_PAGE;
            for(int i=0;i<IND_VERTICES_PER_PAGE ;
                //              && p_vid<num_nodes() ;
                i++,
                v_vid = v_pg*IND_VERTICES_PER_PAGE + i,
                p_vid = p_pg*IND_VERTICES_PER_PAGE + i)
            {
              if(v_vid < prev_num_nodes)
              {//normal copy
                in_delta_vt[p_vid] = orig_delta_vt[i];
              }
              else
              {//fill the rest with null vts
                in_delta_vt[p_vid] =
                  delta_vt_elem(newest_level+1, in_nullcont_idx, 0, 0);
              }
            }
          } //v_pg<prev_pages
          else //new pages. indirections are set to zero
          {
            uint64_t v_vid = v_pg*IND_VERTICES_PER_PAGE;
            uint64_t p_vid = p_pg*IND_VERTICES_PER_PAGE;
            for(int i=0;i<IND_VERTICES_PER_PAGE ;
                i++,
                v_vid = v_pg*IND_VERTICES_PER_PAGE + i,
                p_vid = p_pg*IND_VERTICES_PER_PAGE + i)
            { //fill with null vts
              ASSERT(v_vid ,>=, prev_num_nodes);
              in_delta_vt[p_vid] =
                delta_vt_elem(newest_level+1, in_nullcont_idx, 0, 0);
            }
          }
          //vt cow done. not re-point
          in_indirection[v_pg] =  in_delta_vt + p_pg*IND_VERTICES_PER_PAGE;
          p_pg++;
        }

        //put empty pages
        if(empty_in_pages.size()!=0)
        {
          //create the empty page shared instance
          uint64_t p_vid = p_pg*IND_VERTICES_PER_PAGE;
          for(int i=0;i<IND_VERTICES_PER_PAGE ;
              i++,
              p_vid = p_pg*IND_VERTICES_PER_PAGE + i)
          { //fill with null vts
            in_delta_vt[p_vid] =
              delta_vt_elem(newest_level+1, in_nullcont_idx, 0, 0);
          }

          for(uint64_t v_pg: empty_in_pages)
          { //all to the shared empty page
            ASSERT(v_pg*IND_VERTICES_PER_PAGE ,>=, prev_num_nodes);
            in_indirection[v_pg] = in_delta_vt + p_pg*IND_VERTICES_PER_PAGE;
          }
          p_pg++;
        }
        ASSERT((p_pg*IND_VERTICES_PER_PAGE) ,==, in_delta_size);
      }

/******STEP 4. put the new edges*****/

      {
        uint64_t edge_index=0;
        for(auto elist : out_edges)
        {
          uint64_t start_edge_index=edge_index;
          NodeID_ src = elist.first;
          for(DestID_ dest : elist.second)
          {
            out_et[edge_index++] = dest;
          }
          //put the continuation record
          delta_vt_elem* src_vt =
            out_indirection[IND_V_PG(src)]
            + (IND_V_OFFSET(src));

          *(delta_vt_elem*)(out_et+edge_index) =
            *src_vt;

          //put null cont instead
          //the existing vt entry should already have length zero
          //when doing the cow
          if(src_vt->length() == 0)
          {
            *(delta_vt_elem*)(out_et+edge_index) = CONT_NULL;
          }


          edge_index+=2; //16 bytes

          uint32_t old_degree=src_vt->degree();

          //now replace the vt entry to the new edgetable
          *src_vt =
            delta_vt_elem(newest_level+1, start_edge_index,
                elist.second.size(), old_degree+elist.second.size());

          ASSERT(elist.second.size()+2 ,==,
              edge_index - start_edge_index);
        }

        ASSERT(edge_index, ==,out_nullcont_idx);
        *(delta_vt_elem*)(out_et + out_nullcont_idx) =
          CONT_NULL;
      }




      //now the in edges
      {
        uint64_t edge_index=0;
        for(auto elist : in_edges)
        {
          uint64_t start_edge_index=edge_index;
          NodeID_ src = elist.first;
          for(DestID_ dest : elist.second)
          {
            in_et[edge_index++] = dest;
          }
          //put the continuation record
          delta_vt_elem* src_vt =
            in_indirection[IND_V_PG(src)]
            + (IND_V_OFFSET(src));

          *(delta_vt_elem*)(in_et+edge_index) =
            *src_vt;

          //put null cont instead
          //the existing vt entry should already have length zero
          //when doing the cow
          if(src_vt->length() == 0)
          {
            *(delta_vt_elem*)(in_et+edge_index) = CONT_NULL;
          }


          edge_index+=2; //16 bytes

          uint32_t old_degree=src_vt->degree();

          //now replace the vt entry to the new edgetable
          *src_vt =
            delta_vt_elem(newest_level+1, start_edge_index,
                elist.second.size(), old_degree+elist.second.size());

          ASSERT(elist.second.size()+2 ,==,
              edge_index - start_edge_index);
        }

        ASSERT(edge_index, ==,in_nullcont_idx);
        *(delta_vt_elem*)(in_et + in_nullcont_idx) =
          CONT_NULL;
      }

/****STEP 5. add the delta to the indirection list************/
      out_indirections_.push_back(out_indirection);
      in_indirections_.push_back(in_indirection);

      uint64_t prev_edges = num_edges(newest_level);

      DeltaGraph<NodeID_> dg
        (next_num_nodes,
         prev_edges+el.size(),
         out_delta_vt,  in_delta_vt,
         out_et, in_et,
         out_delta_size, in_delta_size,
         out_et_size, in_et_size
         );
      deltas_.push_back(dg);


      newest_level++;
    }

    bool directed() const
    {
      //all levels should have the identical 'directed' property, right?
      return deltas_[0].directed();
    }

    int64_t out_degree(NodeID_ v, int32_t level=-1) const
    {
      //-1 means 'use the most recent one'
      if (level==-1)
      {
        level=newest_level;
      }
      return (*(out_indirections_[level][IND_V_PG(v)]+IND_V_OFFSET(v))).degree();
//      return deltas_[level].out_degree(v);
    }

    int64_t in_degree(NodeID_ v, int32_t level=-1) const
    {
      //-1 means 'use the most recent one'
      if (level==-1)
      {
        level=newest_level;
      }

      return (*(in_indirections_[level][IND_V_PG(v)]+IND_V_OFFSET(v))).degree();
//      return deltas_[level].in_degree(v);
    }

    //multi-level iterator
    ///it handles the connection between the snapshot, pages
    class ml_iterator
    {
      //uint32_t length_this_level;
      //curr_level means where the iterator is physically on,
      //not the level of your query.
      uint32_t curr_level;
      DestID_* curr_edge;

      typename DeltaGraph<NodeID_>::Neighborhood delta_neigh;

      bool is_out_iter;
      const Snapshots<NodeID_>* snp_;


      public:
      ml_iterator(const Snapshots<NodeID_>* snp):
        snp_(snp){}
      //{snp_=snp;}


      DestID_ operator*() const
      {
        return *curr_edge;
      }



      //in_begin and out_begin : get the curr_edge, curr_level
      //and (base_neigh/delta_neigh) to the
      //beginning of the adjacency list.
      //
      //The codes for these methods are dirty,
      //because I tried to keep a different structure
      //for base graph and delta graph.
      //If I can come up with a elegant way to put
      //abstraction on top of base and delta graphs,
      //I will fix it.
      //
      //maybe i should change them into constructors....

      void in_begin(NodeID_ v, int32_t level=-1)
      {
        is_out_iter = 0;
        int32_t queried_level;
        if(level==-1) queried_level = snp_->newest_level;
        else          queried_level = level;

/****STEP 1. get a vertex entry (watch out for level 0)***/

        delta_vt_elem vertex = *((snp_->in_indirections_[queried_level][IND_V_PG(v)]) + (IND_V_OFFSET(v)));
/****STEP 2. find the starting edge edges (watch out for level 0)*********/
        uint32_t edge_level = vertex.level();
        delta_neigh = snp_->deltas_[edge_level].in_delta_neigh(vertex);
        curr_edge = delta_neigh.begin();
        curr_level = edge_level;
        __builtin_prefetch(curr_edge);

      }


      void out_begin(NodeID_ v, int32_t level=-1)
      {
        is_out_iter = 1;
        int32_t queried_level;
        if(level==-1) queried_level = snp_->newest_level;
        else          queried_level = level;

/****STEP 1. get a vertex entry (watch out for level 0)***/

        delta_vt_elem vertex = *((snp_->out_indirections_[queried_level][IND_V_PG(v)]) + (IND_V_OFFSET(v)));
/****STEP 2. find the starting edge edges (watch out for level 0)*********/
        uint32_t edge_level = vertex.level();
        delta_neigh = snp_->deltas_[edge_level].out_delta_neigh(vertex);
        curr_edge = delta_neigh.begin();
        curr_level = edge_level;
        __builtin_prefetch(curr_edge);

      }

      bool is_end()
      {
        if(
         ((curr_level == 0) && (curr_edge == delta_neigh.end()))
          ||
          ((curr_level!= 0) && (curr_edge == delta_neigh.end())
           && (*curr_edge == CONT_NULL)   )
          )
        {
          return 1;
        }
        else return 0;
      }

      //++ the iterator, and follow the cont. record.
      //return whether it is at the end? (need it?)
      //is there any elegant way of doing this?
      //
      bool next()
      {
        curr_edge++;
        //follow cont.record. no need to handle curr_level==0 case
        if((curr_level != 0) && (curr_edge == delta_neigh.end()))
        {
          //cont.record is in the form of delta_vt_elem. so we can follow.
          delta_vt_elem* cont = (delta_vt_elem*)curr_edge;

          if(!cont->cont_record_null()) //traverse more levels
          {
            curr_level = cont->level();
            if(is_out_iter)
            {
              delta_neigh = snp_->deltas_[curr_level].
                out_delta_neigh(*cont);
            } else
            {
              delta_neigh = snp_->deltas_[curr_level].
                in_delta_neigh(*cont);
            }
            curr_edge = delta_neigh.begin();
          }
        }
        else
        {
          __builtin_prefetch(curr_edge+32);
        }
        return is_end();
      }//next()
    };//ml_iterator


    ml_iterator iter() const
    {
      return ml_iterator(this);
    }

    virtual void PrintStats() const {
      if(newest_level == -1) return;
      std::cout << "Graph has " << num_nodes() << " nodes and "
        << num_edges() << " ";
      if (!directed())
        std::cout << "un";
      std::cout << "directed edges for avg degree: ";
      std::cout << num_edges()/num_nodes() << std::endl;
    }


  private:
    uint64_t count_edges(EdgeList &el, uint64_t &next_num_nodes,
        std::set<uint64_t> &new_out_pages ,
        std::set<uint64_t> &new_in_pages ,
      std::map<NodeID_,std::list<DestID_>> &out_edges,
      std::map<NodeID_,std::list<DestID_>> &in_edges)
    {
      int next_num_edges=0;
      for(auto e : el)
      {
        next_num_edges++;
        new_out_pages.insert(IND_V_PG(e.u));
        new_in_pages.insert(IND_V_PG(e.v));

        //new_out_vertices.insert(e.u);
        //new_in_vertices.insert(e.v);

        out_edges[e.u].push_back(e.v);
        in_edges[e.v].push_back(e.u);


        if(e.u >= next_num_nodes)
        {
          next_num_nodes = e.u+1;
        }
        if(e.v >= next_num_nodes)
        {
          next_num_nodes = e.v+1;
        }
      }

      return next_num_edges;
    }


    uint64_t count_empty_pages(uint64_t prev_pages, uint64_t next_pages,
        std::set<uint64_t> new_out_pages ,
        std::set<uint64_t> new_in_pages ,
      std::set<uint64_t>& empty_out_pages,
      std::set<uint64_t>& empty_in_pages)
    {
      for(uint64_t page= prev_pages;page < next_pages; page++)
      {
        if(!new_out_pages.count(page))
        {
          empty_out_pages.insert(page);
        }
        if(!new_in_pages.count(page))
        {
          empty_in_pages.insert(page);
        }
      }
      return next_pages;
    }


/*****STEP 2. alloc the structures*******************/
      //input : new_io_pages  empty_io_pages dbpath next_num_edges num_cont_record
      //output: io_delta_size(dbg)  io_delta_vt  io_et io_nullcont_idx
    void alloc_tables(
        //inputs
      std::set<uint64_t> new_out_pages,
      std::set<uint64_t> new_in_pages,
      std::set<uint64_t> empty_out_pages,
      std::set<uint64_t> empty_in_pages ,
      std::map<NodeID_,std::list<DestID_>> out_edges,
      std::map<NodeID_,std::list<DestID_>> in_edges ,
      uint64_t next_num_edges,
      std::string dbpath,
      bool is_base,//to determine whether we need cont.records
      //outputs
      delta_vt_elem* *out_delta_vt,
      delta_vt_elem* *in_delta_vt ,
      DestID_* *out_et,
      DestID_* *in_et ,
      uint64_t* p_out_delta_size,
      uint64_t* p_in_delta_size,
      uint64_t* p_out_et_size,
      uint64_t*  p_in_et_size,
      uint64_t* in_nullcont_idx ,
      uint64_t* out_nullcont_idx
        )
    {
      uint64_t out_delta_size= IND_VERTICES_PER_PAGE*new_out_pages.size();
      if(empty_out_pages.size()!=0)
      {//for the shared empty page
        out_delta_size+=IND_VERTICES_PER_PAGE;
      }
      uint64_t in_delta_size= IND_VERTICES_PER_PAGE*new_in_pages.size();
      if(empty_in_pages.size()!=0)
      {
        in_delta_size+=IND_VERTICES_PER_PAGE;
      }

      *out_delta_vt =
        (delta_vt_elem*)mmap_alloc("out_V_", dbpath,newest_level+1,
          sizeof(delta_vt_elem) * out_delta_size);
      //std::cout<<"outdt alloced at "<<*out_delta_vt<<" with size "<<out_delta_size<<" (sizeof(delta_vt_elem)*out_delta_size "<<std::hex<<sizeof(delta_vt_elem) * out_delta_size<<" = "<<*out_delta_vt + out_delta_size<<std::endl;

      *in_delta_vt =
        (delta_vt_elem*)mmap_alloc("in_V_", dbpath, newest_level+1,
          sizeof(delta_vt_elem) * in_delta_size);

//      std::cout<<"indt alloced at "<<*in_delta_vt<<" with size "<<in_delta_size<<" (sizeof(delta_vt_elem)*in_delta_size "<<std::hex<<sizeof(delta_vt_elem) * in_delta_size<<" = "<<*in_delta_vt + in_delta_size<<std::endl;


      //edges plus cont.records( plus one nullcont)

      uint64_t out_et_size;
      uint64_t in_et_size ;

      if(is_base)
      {
          out_et_size =(next_num_edges + CONT_SIZE );
          in_et_size = (next_num_edges + CONT_SIZE );

          //nullcont will be stored here
          *in_nullcont_idx = next_num_edges;
          *out_nullcont_idx = next_num_edges;
      }
      else//delta
      {

        out_et_size =
          (next_num_edges + (1+out_edges.size()) * CONT_SIZE);
        in_et_size =
          (next_num_edges + (1+in_edges.size()) * CONT_SIZE);

        //nullcont will be stored here
        *in_nullcont_idx = next_num_edges +
          in_edges.size() * CONT_SIZE;
        *out_nullcont_idx = next_num_edges +
          out_edges.size() * CONT_SIZE;
      }

      *p_out_delta_size= out_delta_size;
      *p_in_delta_size=  in_delta_size;
      *p_out_et_size= out_et_size;
      *p_in_et_size=  in_et_size;


      *out_et =
        (DestID_*)mmap_alloc("out_E_", dbpath, newest_level+1,
            sizeof(DestID_) * out_et_size);

      //      std::cout<<"in edge alloced at "<<*out_et<<" with size "<<out_et_size<<" (sizeof(DestID_)*next_num_ed "<<std::hex<<sizeof(DestID_) * next_num_edges<<" = "<<*out_et + next_num_edges<<std::endl;

      *in_et =
        (DestID_*)mmap_alloc("in_E_", dbpath, newest_level+1 ,
            sizeof(DestID_) * in_et_size);

      //      std::cout<<"in edge alloced at "<<*in_et<<" with size "<<in_et_size<<" (sizeof(DestID_)*next_num_ed "<<std::hex<<sizeof(DestID_) * next_num_edges<<" = "<<*in_et + next_num_edges<<std::endl;
    }
    std::vector<int64_t>  out_neigh(int64_t u) const
    {
        ml_iterator miter = this->iter();
        std::vector<int64_t> out_neigh;
        for(miter.out_begin(u);!miter.is_end();miter.next())
        {
            Add a comment to this line
                out_neigh.push_back(*miter);

        }

        return out_neigh;

    }


    std::vector<int64_t>  in_neigh(int64_t u) const
    {
        ml_iterator miter = this->iter();
        std::vector<int64_t> in_neigh;
        for(miter.in_begin(u);!miter.is_end();miter.next())
        {
            in_neigh.push_back(*miter);

        }

        return in_neigh;

    }














};





#endif //SNAPSHOT_H_

