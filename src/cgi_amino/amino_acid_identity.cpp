/**
 * @file    amino_acid_identity.cpp
 * @ingroup src
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include <iostream>
#include <ctime>
#include <chrono>
#include <functional>

//Own includes
#include "cgi_amino/include/aai_parameters.hpp"
#include "cgi_amino/include/parseCmdArgs.hpp"
#include "cgi_amino/include/computeIdentity.hpp"
#include "map/include/base_types.hpp"

//External includes
#include "common/argvparser.hpp"

int main(int argc, char** argv)
{
  CommandLineProcessing::ArgvParser cmd;

  //Setup command line options
  aa::initCmdParser(cmd);

  //Parse command line arguements   
  aa::Parameters parameters;       
  aa::parseandSave(argc, argv, cmd, parameters);   

  auto t0 = skch::Time::now();

  aa::Compute aai_solver(parameters);

  std::chrono::duration<double> timeAAIcompute = skch::Time::now() - t0;
  std::cerr << "INFO, aai::main, Time spent computing AAI : " << timeAAIcompute.count() << " sec" << std::endl;
}
