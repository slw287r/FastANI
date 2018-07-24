/**
 * @file    core_genome_identity.cpp
 * @ingroup src
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include <iostream>
#include <ctime>
#include <chrono>
#include <functional>

//Own includes
#include "map/include/map_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/parseCmdArgs.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/computeMap.hpp"
#include "map/include/commonFunc.hpp"
#include "cgi/include/computeCoreIdentity.hpp" 

//External includes
#include "common/argvparser.hpp"

int main(int argc, char** argv)
{
  CommandLineProcessing::ArgvParser cmd;
  using namespace std::placeholders;  // for _1, _2, _3...

  //Setup command line options
  skch::initCmdParser(cmd);

  //Parse command line arguements   
  skch::Parameters parameters;        //sketching and mapping parameters

  skch::parseandSave(argc, argv, cmd, parameters);   

  //Redirect Mashmap's mapping output to null fs, using file name for CGI output
  std::string fileName = parameters.outFileName;

#ifdef DEBUG
  parameters.outFileName = parameters.outFileName + ".map";
#else
  parameters.outFileName = "/dev/null";
#endif

  auto t0 = skch::Time::now();

  //Build the sketch for reference
  skch::Sketch referSketch(parameters);

  std::chrono::duration<double> timeRefSketch = skch::Time::now() - t0;
  std::cerr << "INFO, skch::main, Time spent sketching the reference : " << timeRefSketch.count() << " sec" << std::endl;

  //Initialize the files to delete the existing content
  {
#ifdef DEBUG
    std::ofstream outstrm2(fileName + ".map.1way");
    std::ofstream outstrm3(fileName + ".map.2way");
#endif
    if(parameters.visualize)
      std::ofstream outstrm4(fileName + ".visual");
  }

  //Final output vector of ANI computation
  std::vector<cgi::CGI_Results> finalResults;

  //Loop over query genomes
  for(uint64_t queryno = 0; queryno < parameters.querySequences.size(); queryno++)
  {
    t0 = skch::Time::now();

    skch::MappingResultsVector_t mapResults;
    uint64_t totalQueryFragments = 0;

    auto fn = std::bind(skch::Map::insertL2ResultsToVec, std::ref(mapResults), _1);
    skch::Map mapper = skch::Map(parameters, referSketch, totalQueryFragments, queryno, fn);

    std::chrono::duration<double> timeMapQuery = skch::Time::now() - t0;
    std::cerr << "INFO, skch::main, Time spent mapping fragments in query #" << queryno + 1 <<  " : " << timeMapQuery.count() << " sec" << std::endl;

    t0 = skch::Time::now();

    cgi::computeCGI(parameters, mapResults, mapper, referSketch, totalQueryFragments, queryno, fileName, finalResults);

    std::chrono::duration<double> timeCGI = skch::Time::now() - t0;
    std::cerr << "INFO, skch::main, Time spent post mapping : " << timeCGI.count() << " sec" << std::endl;
  }

  cgi::outputCGI (parameters, finalResults, fileName);

  if (parameters.matrixOutput)
    cgi::outputPhylip (parameters, finalResults, fileName);
}
