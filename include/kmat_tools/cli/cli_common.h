#pragma once

#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include <bcli/bcli.hpp>

namespace kmat
{

enum class COMMAND
{
  CONVERT,
  DIFF,
  FASTA,
  FAFMT,
  FILTER,
  KTFILTER,
  MERGE,
  REVERSE,
  SELECT,
  UNITIG,
  UNKNOWN
};

inline COMMAND str_to_cmd(const std::string& s)
{
  if (s == "convert")
    return COMMAND::CONVERT;
  else if (s == "diff")
    return COMMAND::DIFF;
  else if (s == "fasta")
    return COMMAND::FASTA;
  else if (s == "fafmt")
    return COMMAND::FAFMT;
  else if (s == "filter")
    return COMMAND::FILTER;
  else if (s == "ktfilter")
    return COMMAND::KTFILTER;
  else if (s == "merge")
    return COMMAND::MERGE;
  else if (s == "reverse")
    return COMMAND::SELECT;
  else if (s == "unitig")
    return COMMAND::UNITIG;
  else 
    return COMMAND::UNKNOWN;
}

inline std::string cmd_to_str(COMMAND cmd)
{
  if (cmd == COMMAND::CONVERT)
    return "convert";
  else if (cmd == COMMAND::DIFF)
    return "diff";
  else if (cmd == COMMAND::FASTA)
    return "fasta";
  else if (cmd == COMMAND::FAFMT)
    return "fafmt";
  else if (cmd == COMMAND::FILTER)
    return "filter";
  else if (cmd == COMMAND::KTFILTER)
    return "ktfilter";
  else if (cmd == COMMAND::MERGE)
    return "merge";
  else if (cmd == COMMAND::SELECT)
    return "select";
  else if (cmd == COMMAND::UNITIG)
    return "unitig";
  else // (cmd == COMMAND::UNKNOWN)
    return "unknown";
}

struct kmat_options { };
using kmat_opt_t = std::shared_ptr<struct kmat_options>;

using cli_t = std::shared_ptr<bc::Parser<1>>;

};