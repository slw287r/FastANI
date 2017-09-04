/**
 * @file    aai_types.hpp
 * @brief   specific type definitions for AAI computation
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef AAI_BASE_TYPES_HPP 
#define AAI_BASE_TYPES_HPP

namespace aa
{

  typedef skch::hash_t hash_t;    //hash type
  typedef skch::offset_t offset_t;       //position within sequence
  typedef skch::seqno_t seqno_t;        //sequence counter in file

  //Information about each sketch element value
  struct sketchElementInfo
  {
    hash_t hash;
    seqno_t geneId;

    static bool lessByHash(const sketchElementInfo &x, const sketchElementInfo &y) {
      return x.hash < y.hash;
    }
  };

  //Pair to denote which two genes share a sketch element
  struct geneMatchPairInfo
  {
    seqno_t geneIdR;                      //gene id from reference genome
    seqno_t geneIdQ;                      //gene id from query genome

    static bool lessByRefGene(const geneMatchPairInfo &x, const geneMatchPairInfo &y) {
      return x.geneIdR < y.geneIdR;
    }

    static bool lessByQryGene(const geneMatchPairInfo &x, const geneMatchPairInfo &y) {
      return x.geneIdQ < y.geneIdQ;
    }

    static bool lessByRefAndQryGene(const geneMatchPairInfo &x, const geneMatchPairInfo &y) {
      return std::tie(x.geneIdR, x.geneIdQ) < std::tie(y.geneIdR, y.geneIdQ);
    }

    static bool lessByQryAndRefGene(const geneMatchPairInfo &x, const geneMatchPairInfo &y) {
      return std::tie(x.geneIdQ, x.geneIdR) < std::tie(y.geneIdQ, y.geneIdR);
    }
  };

  struct geneBestMatchInfo
  {
    seqno_t geneIdR;
    seqno_t geneIdQ;
    offset_t countSharedSketchVals;

    bool operator ==(const geneBestMatchInfo &rhs)
    {
      return std::tie(this->geneIdR, this->geneIdQ, this->countSharedSketchVals) == 
        std::tie(rhs.geneIdR, rhs.geneIdQ, rhs.countSharedSketchVals);
    }

    bool operator <(const geneBestMatchInfo &rhs) {
      return std::tie(this->geneIdR, this->geneIdQ, this->countSharedSketchVals) < 
        std::tie(rhs.geneIdR, rhs.geneIdQ, rhs.countSharedSketchVals);
    }
  };
}

#endif
