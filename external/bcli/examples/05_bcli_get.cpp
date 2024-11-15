#include<bcli/bcli.hpp>

using namespace bc;

class UserAddr
{
public:
  UserAddr(const std::string& user, const std::string& addr)
    : user(user), addr(addr)
  {}
  void show() {std::cerr << "user:" << user << ", addr:" << addr << std::endl;}
private:
  std::string user;
  std::string addr;
};

class Simple
{
public:
  Simple(const std::string& v) : v(v) {}
  std::string v;
};

template<>
UserAddr param::Param::as<UserAddr>()
{
  std::vector<std::string> s = utils::split(m_str_value, '@');
  return UserAddr(s[0], s[1]);
}

int main(int argc, char* argv[])
{
  Parser<0> cli("ex05",
                "ex05 desc, try: ./05_bcli_settter -s 10",
                "v0.0.1",
                "J. Doe");

  cli.add_param("-s", "string")->def("txt");
  cli.add_param("-i", "int")->def("32");
  cli.add_param("-f", "float")->def("18.6");
  cli.add_param("-a", "user@addr")->def("user@localhost");

  cli.add_common();

  BCLI_PARSE(cli, argc, argv)

  std::string s = cli.getp("s")->as<std::string>();
  int i = cli.getp("i")->as<int>();
  int f = cli.getp("f")->as<float>();

  UserAddr A = cli.getp("a")->as<UserAddr>();
  Simple   S = cli.getp("s")->as<Simple>();
  
  std::cerr << "s " << s << std::endl;
  std::cerr << "i " << i << std::endl;
  std::cerr << "f " << f << std::endl;
  std::cerr << "Simple " << S.v << std::endl;
  A.show();

}

// ex05 v0.0.1
// 
// DESCRIPTION
//   ex05 desc, try: ./05_bcli_settter -s 10
// 
// USAGE
//   ex05 [-s <?>] [-i <?>] [-f <?>] [-a <?>] [-h/--help] [-v/--verbose] [-d/--debug] [--version] 
// 
// OPTIONS
//   [global] - global parameters
//     -s  - string {txt}
//     -i  - int {32}
//     -f  - float {18.6}
//     -a  - user@addr {user@localhost}
// 
//   [common]
//     -h --help    - Show this message and exit. [⚑]
//     -v --verbose - Verbose mode. [⚑]
//     -d --debug   - Debug mode. [⚑]
//        --version - Show version and exit. [⚑]