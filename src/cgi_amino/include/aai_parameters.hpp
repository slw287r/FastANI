/**
 * @file    aai_parameters.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef AAI_CONFIG_HPP 
#define AAI_CONFIG_HPP

#include <vector>

namespace aa
{
  /**
   * @brief   configuration parameters for computing AAI
   *          expected to be initialized using command line arguments
   */
  struct Parameters
  {
    int kmerSize;                                     //kmer size for sketching
    int minGeneLength;                                //minimum gene size to use for AAI computation
    int sketchSize;                                   //sketch size per gene
    int minFragments;                                 //minimum gene mappings for trusting AAI value
    int alphabetSize;                                 //alphabet size
    std::vector<std::string> refSequences;            //reference sequence(s)
    std::vector<std::string> querySequences;          //query sequence(s)
    std::string outFileName;                          //output file name
  };

  /**
   * @brief     Internal figures not exposed at the command line interface
   */
  namespace fixed
  {
    float minimumIdentity = 30.0;                     //minimum identity between bi-directional mappings
  }
}

#endif
