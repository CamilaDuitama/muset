#define _CFRCAT_APP_
#include <cfrcat/cfrcat.hpp>
#include <cfrcat/options.hpp>

int main(int argc, char* argv[])
{
  cfr::cfrcat_params params = cfr::parse(argc, argv);
  cfr::process(params);
  return 0;
}