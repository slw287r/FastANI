/**
 * @file    computeIdentity.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMPUTE_IDENTITY_HPP 
#define COMPUTE_IDENTITY_HPP

#include <vector>
#include <algorithm>

//Own includes
#include "cgi_amino/include/sketch.hpp"


namespace aa
{
  class Compute
  {
    //algorithm parameters
    const aa::Parameters &param;

    public:

    /**
     * @brief   constructor
     *          also calls routine for distance computation
     */
    Compute(const aa::Parameters &p) 
      :
        param(p)
        {
          this->computeAAI();
        }


    private:

    /**
     * @brief     main routine to compute AAI between all given genomes
     */
    void computeAAI()
    {
      std::ofstream outstrm(param.outFileName);

      //If query and reference genomes are same,
      //avoid renundant combinations

      if(param.refSequences != param.querySequences)
      {
      for(const auto &rFileName : param.refSequences)
        for(const auto &qFileName : param.querySequences)
        {
#ifdef DEBUG
          std::cerr << "INFO, aa::Compute::computeAAI, computing AAI for reference genome =" << rFileName
            << " and query genome =" << qFileName << std::endl;
#endif
          this->computeSinglePair(rFileName, qFileName, outstrm);
        }
      }
      else
      {
        int loopSize = param.refSequences.size();

        for(int i = 0; i < loopSize; i++)
          for(int j = i; j < loopSize; j++)
          {
#ifdef DEBUG
          std::cerr << "INFO, aa::Compute::computeAAI, computing AAI for reference genome =" << rFileName
            << " and query genome =" << qFileName << std::endl;
#endif
          this->computeSinglePair(param.refSequences[i], param.querySequences[j], outstrm);
          }
      }
    }

    /**
     * @brief                       compute AAI between single pair of genomes
     * @param[in]   rFileName       reference genome file 
     * @param[in]   qFileName       query genome file
     * @param[in]   outstrm         output stream object
     */
    void computeSinglePair(const std::string &rFileName, const std::string &qFileName, std::ofstream &outstrm)
    {
      //built sketch index for given genomes
      Sketch refSkch(this->param, rFileName);
      Sketch qrySkch(this->param, qFileName);
      
      //compute pairs of genes that share any sketch element 
      //number of same pairs indicate the count of shared elements
      std::vector< geneMatchPairInfo > geneMatchPairs; 
      this->findGenePairs(refSkch, qrySkch, geneMatchPairs); 

      //compute reciprocal best gene hits
      std::vector< geneBestMatchInfo > reciprocalGeneMatches;
      this->computeGeneMappings(geneMatchPairs, reciprocalGeneMatches);

      this->computeAndReportAAI(reciprocalGeneMatches, rFileName, qFileName, outstrm);
    }

    /**
     * @brief                         compute pairs of genes that share any sketch element           
     * @param[in]   r                 reference genome sketch index
     * @param[in]   q                 query genome sketch index
     * @param[out]  geneMatchPairs    output stream object
     */
    template <typename VEC>
      void findGenePairs(Sketch &r, Sketch &q, VEC &geneMatchPairs)
      {
        std::sort(r.sketchIndex.begin(), r.sketchIndex.end(), sketchElementInfo::lessByHash);
        std::sort(q.sketchIndex.begin(), q.sketchIndex.end(), sketchElementInfo::lessByHash);

        auto it_q = q.sketchIndex.begin(); 

        for(auto it_r = r.sketchIndex.begin(); it_r != r.sketchIndex.end();)
        {
          auto it_r_rangeEnd = std::upper_bound(it_r, r.sketchIndex.end(), *it_r, sketchElementInfo::lessByHash);

          it_q = std::lower_bound(it_q, q.sketchIndex.end(), *it_r, sketchElementInfo::lessByHash); 
          auto it_q_rangeEnd =  std::upper_bound(it_q, q.sketchIndex.end(), *it_r, sketchElementInfo::lessByHash);

          //If reference sketch value exists in query, then [it_q, it_q_rangeEnd) represents the range in query genome 
          //with same value

          for(auto i = it_q; i != it_q_rangeEnd; i++)
            for(auto j = it_r; j != it_r_rangeEnd; j++)
            {
              //Push reference gene id, query gene id into vector
              geneMatchPairs.emplace_back( geneMatchPairInfo{j->geneId, i->geneId} );
            }

          it_r = it_r_rangeEnd;
        }
      }

    /**
     * @brief                               compute bi-directional gene mappings and count 
     *                                      of shared sketch elements between them
     * @param[in]   geneMatchPairs          sketch elements sharing between genes
     * @param[out]  reciprocalGeneMatches   gene mappings
     */
    template <typename VEC1, typename VEC2>
      void computeGeneMappings(VEC1 &geneMatchPairs, VEC2 &reciprocalGeneMatches)
      {
        std::vector< geneBestMatchInfo > bestGeneMatchesForRef;

        {
          std::sort(geneMatchPairs.begin(), geneMatchPairs.end(), geneMatchPairInfo::lessByRefAndQryGene);
          
          for(auto it_r = geneMatchPairs.begin(); it_r != geneMatchPairs.end(); )
          {
            //Range with same reference gene id
            auto it_r_end = std::upper_bound(it_r, geneMatchPairs.end(), *it_r, geneMatchPairInfo::lessByRefGene);

            seqno_t bestQueryGeneId;
            int maxFrequency = 0;

            //Now find maximum occuring query geneid in this range
            for(auto it_q = it_r; it_q != it_r_end;)
            {
              auto it_q_end = std::upper_bound(it_q, it_r_end, *it_q, geneMatchPairInfo::lessByQryGene);

              if(std::distance(it_q, it_q_end) > maxFrequency)
              {
                bestQueryGeneId = it_q->geneIdQ;
                maxFrequency = std::distance(it_q, it_q_end);
              }

              it_q = it_q_end;
            }

            bestGeneMatchesForRef.emplace_back( geneBestMatchInfo{it_r->geneIdR, bestQueryGeneId, maxFrequency}  ); 

            //Advance the iterator
            it_r = it_r_end;
          }
        }

        std::vector< geneBestMatchInfo > bestGeneMatchesForQry;

        {
          std::sort(geneMatchPairs.begin(), geneMatchPairs.end(), geneMatchPairInfo::lessByQryAndRefGene);
          
          for(auto it_q = geneMatchPairs.begin(); it_q != geneMatchPairs.end(); )
          {
            //Range with same reference gene id
            auto it_q_end = std::upper_bound(it_q, geneMatchPairs.end(), *it_q, geneMatchPairInfo::lessByQryGene);

            seqno_t bestRefGeneId;
            int maxFrequency = 0;

            //Now find maximum occuring query geneid in this range
            for(auto it_r = it_q; it_r != it_q_end;)
            {
              auto it_r_end = std::upper_bound(it_r, it_q_end, *it_r, geneMatchPairInfo::lessByRefGene);

              if(std::distance(it_r, it_r_end) > maxFrequency)
              {
                bestRefGeneId = it_r->geneIdR;
                maxFrequency = std::distance(it_r, it_r_end);
              }

              it_r = it_r_end;
            }

            bestGeneMatchesForQry.emplace_back( geneBestMatchInfo{bestRefGeneId, it_q->geneIdQ, maxFrequency}  ); 

            //Advance the iterator
            it_q = it_q_end;
          }

        }

        std::vector< geneBestMatchInfo > bestGeneMatchesBoth;

        {
          bestGeneMatchesBoth.insert(bestGeneMatchesBoth.end(), bestGeneMatchesForRef.begin(), bestGeneMatchesForRef.end());
          bestGeneMatchesBoth.insert(bestGeneMatchesBoth.end(), bestGeneMatchesForQry.begin(), bestGeneMatchesForQry.end());

          std::sort(bestGeneMatchesBoth.begin(), bestGeneMatchesBoth.end());

          for(auto it = bestGeneMatchesBoth.begin(); it != bestGeneMatchesBoth.end(); it++)
          {
            if( std::next(it) == bestGeneMatchesBoth.end() )
              break;

            auto it2 = std::next(it);

            if( *it == *it2 )
              reciprocalGeneMatches.emplace_back( geneBestMatchInfo{ it->geneIdR, it->geneIdQ, it->countSharedSketchVals } );
          }
        }
      }

    /**
     * @brief                               compute and report AAI using reciprocal mappings
     * @param[in]   reciprocalGeneMatches   gene mappings
     * @param[in]   rFileName               reference genome file 
     * @param[in]   qFileName               query genome file
     * @param[in]   outstrm                 output stream object
     */
    template <typename VEC1>
      void computeAndReportAAI(VEC1 &reciprocalGeneMatches, 
          const std::string &rFileName, const std::string &qFileName,
          std::ofstream &outstrm)
      {
        int countOfMappings = 0;
        float sumIdentity = 0;

        for(auto &e : reciprocalGeneMatches)
        {
          float jaccard = e.countSharedSketchVals * 1.0 / param.sketchSize;
          float distance_estimate = skch::Stat::j2md(jaccard, param.kmerSize);
          float identity = 100 * (1 - distance_estimate);

          if(identity >= fixed::minimumIdentity) 
          {
            sumIdentity += identity;
            countOfMappings += 1;
          }
        }

        if(countOfMappings >= param.minFragments)
        {
          //Report AAI
          outstrm << rFileName 
            << " " << qFileName
            << " " << (countOfMappings > 0 ? sumIdentity/countOfMappings : 0.0)
            << " " << countOfMappings 
#ifdef DEBUG
            << std::endl;
#else
            << "\n";
#endif
        }
      }
  };
}

#endif
