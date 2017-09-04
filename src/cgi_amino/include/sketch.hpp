/**
 * @file    sketch.hpp
 * @brief   Routines to sketch genomes (protein annotations) 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef BUILD_SKETCH_HPP 
#define BUILD_SKETCH_HPP

#include <zlib.h>  
#include <vector>
#include <algorithm>

//Own includes
#include "cgi_amino/include/aai_types.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/kseq.h"
#include "common/murmur3.h"

KSEQ_INIT(gzFile, gzread)

namespace aa
{
  /**
   * @class     aa::Sketch
   * @brief     computes sketch of the given genome
   */
  class Sketch
  {
    //private members
    
    //seed for murmerhash
    const int seed = 42;

    //algorithm parameters
    const aa::Parameters &param;

    const std::string &fileName;

    public:

    typedef std::vector< sketchElementInfo > MI_Type;
    using MIIter_t = MI_Type::const_iterator;

    MI_Type sketchIndex;

    /**
     * @brief   constructor
     *          also builds the sketch index
     */
    Sketch(const aa::Parameters &p, const std::string &file) 
      :
        param(p),
        fileName(file) 
        {
          this->build();
        }

    private:

    void build()
    {
      FILE *file = fopen(fileName.c_str(), "r");
      gzFile fp = gzdopen(fileno(file), "r");
      kseq_t *seq = kseq_init(fp);

      //sequence counter while parsing file
      seqno_t geneCounter = 0;

      //size of sequence
      offset_t len;

      while ((len = kseq_read(seq)) >= 0) 
      {
        //Is the sequence too short?
        if(len < param.kmerSize || len < param.minGeneLength)
        {
          //do nothing
        }
        else
        {
          this->computeSketch(seq, geneCounter);
        }

        geneCounter++;
      }

      kseq_destroy(seq);  
      gzclose(fp); //close the file handler 
      fclose(file);
    }

    template <typename KSEQ>
      void computeSketch(KSEQ kseq, seqno_t geneCounter)
      {
        //length of the sequence
        offset_t len = kseq->seq.l;

        std::vector<hash_t> sketchElmentValues;

        //Parse all kmers in given sequence
        for(offset_t i = 0; i < len - param.kmerSize + 1; i++)
        {
          hash_t sketchVal = skch::CommonFunc::getHash(kseq->seq.s + i, param.kmerSize); 

          //we don't need to consider opposite strand because this is protein sequence

          sketchElmentValues.emplace_back(sketchVal);    
        }

        std::sort(sketchElmentValues.begin(), sketchElmentValues.end());  

        //Filter out the duplicate values
        auto last = std::unique(sketchElmentValues.begin(), sketchElmentValues.end());
        sketchElmentValues.erase(last, sketchElmentValues.end());

        //Only keep minimum 's' values
        if(sketchElmentValues.size() > param.sketchSize)
          sketchElmentValues.resize(param.sketchSize);

        //Put the sketch set into index
        std::for_each(sketchElmentValues.begin(), sketchElmentValues.end(), [&](hash_t &e){
            sketchIndex.emplace_back( sketchElementInfo{e, geneCounter} );
            });
      }

  };
}

#endif
