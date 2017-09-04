/**
 * @file    parseCmdArgs.hpp
 * @brief   Functionality related to command line parsing 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef AAI_PARSE_CMD_HPP 
#define AAI_PARSE_CMD_HPP

#include <iostream>
#include <string>
#include <fstream>

//Own includes
#include "cgi_amino/include/aai_parameters.hpp"
#include "map/include/parseCmdArgs.hpp"

//External includes
#include "common/argvparser.hpp"

namespace aa
{

  /**
   * @brief           Initialize the command line argument parser 
   * @param[out] cmd  command line parser object
   */
  void initCmdParser(CommandLineProcessing::ArgvParser &cmd)
  {
    cmd.setIntroductoryDescription("-----------------\n\
fastAAI is a fast alignment-free implementation for computing AAI between genomes (protein annotations)\n\
-----------------\n\
Example usage: \n\
$ fastAAI -s genome1.faa -q genome2.faa -o output.txt\n\
$ fastAAI --sl genome_list.txt -q genome2.faa -o output.txt");

    cmd.setHelpOption("h", "help", "Print this help page");

    cmd.defineOption("subject", "an input reference file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("subject","s");

    cmd.defineOption("subjectList", "a file containing list of reference genome files, one genome per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("subjectList","sl");

    cmd.defineOption("query", "an input query file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("query","q");

    cmd.defineOption("queryList", "a file containing list of query genome files, one genome per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("queryList","ql");

    cmd.defineOption("kmer", "kmer size <= 7 [default 7]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("kmer","k");

    cmd.defineOption("minFrag", "minimum mappings for trusting AAI [default : 50]", ArgvParser::OptionRequiresValue);

    cmd.defineOption("minGeneLength", "minimum gene length to be used for AAI computation [default : 200]", ArgvParser::OptionRequiresValue);

    cmd.defineOption("sketchSize", "sketch size per gene [default : 30]", ArgvParser::OptionRequiresValue);

    cmd.defineOption("output", "output file name", ArgvParser::OptionRequired | ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("output","o");
  }

  /**
   * @brief                   Print the parsed cmd line options
   * @param[in]  parameters   parameters parsed from command line
   */
  void printCmdOptions(aa::Parameters &parameters)
  {
    std::cerr << ">>>>>>>>>>>>>>>>>>" << std::endl;
    std::cerr << "Reference = " << parameters.refSequences << std::endl;
    std::cerr << "Query = " << parameters.querySequences << std::endl;
    std::cerr << "Kmer size = " << parameters.kmerSize << std::endl;
    std::cerr << "Minimum gene length = " << parameters.minGeneLength << std::endl;
    std::cerr << "Sketch size per gene  = " << parameters.sketchSize << std::endl;
    std::cerr << "Minimum bi-directional mappings  = " << parameters.minFragments << std::endl;
    std::cerr << "AAI output file = " << parameters.outFileName << std::endl;
    std::cerr << ">>>>>>>>>>>>>>>>>>" << std::endl;
  }

  /**
   * @brief                   Parse the cmd line options
   * @param[in]   cmd
   * @param[out]  parameters  sketch parameters are saved here
   */
  void parseandSave(int argc, char** argv, 
      CommandLineProcessing::ArgvParser &cmd, 
      aa::Parameters &parameters)
  {
    int result = cmd.parse(argc, argv);

    //Make sure we get the right command line args
    if (result != ArgvParser::NoParserError)
    {
      std::cerr << cmd.parseErrorDescription(result) << "\n";
      exit(1);
    }
    else if (!cmd.foundOption("subject") && !cmd.foundOption("subjectList"))
    {
      std::cerr << "Provide reference file (s)\n";
      exit(1);
    }
    else if (!cmd.foundOption("query") && !cmd.foundOption("queryList"))
    {
      std::cerr << "Provide query file (s)\n";
      exit(1);
    }

    std::stringstream str;

    //Parse reference files
    if(cmd.foundOption("subject"))
    {
      std::string ref;

      str << cmd.optionValue("subject");
      str >> ref;

      parameters.refSequences.push_back(ref);
    }
    else //list of files
    {
      std::string listFile;

      str << cmd.optionValue("subjectList");
      str >> listFile;

      skch::parseFileList(listFile, parameters.refSequences);
    }

    //Size of reference
    str.clear();

    //Parse query files
    if(cmd.foundOption("query"))
    {
      std::string query;

      str << cmd.optionValue("query");
      str >> query;

      parameters.querySequences.push_back(query);
    }
    else //list of files
    {
      std::string listFile;

      str << cmd.optionValue("queryList");
      str >> listFile;

      skch::parseFileList(listFile, parameters.querySequences);
    }
    
    str.clear();

    parameters.alphabetSize = 20;

    //Parse algorithm parameters
    if(cmd.foundOption("kmer"))
    {
      str << cmd.optionValue("kmer");
      str >> parameters.kmerSize;
      str.clear();
    }
    else
    {
      parameters.kmerSize = 7;
    }

    if(cmd.foundOption("minGeneLength"))
    {
      str << cmd.optionValue("minGeneLength");
      str >> parameters.minGeneLength;
      str.clear();
    }
    else
      parameters.minGeneLength = 200;

    if(cmd.foundOption("minFrag"))
    {
      str << cmd.optionValue("minFrag");
      str >> parameters.minFragments;
      str.clear();
    }
    else
      parameters.minFragments = 50;

    if(cmd.foundOption("sketchSize"))
    {
      str << cmd.optionValue("sketchSize");
      str >> parameters.sketchSize;
      str.clear();
    }
    else
      parameters.sketchSize = 30;


    str << cmd.optionValue("output");
    str >> parameters.outFileName;
    str.clear();

    printCmdOptions(parameters);

    //Check if files are valid
    skch::validateInputFiles(parameters.querySequences, parameters.refSequences);
  }
}

#endif
