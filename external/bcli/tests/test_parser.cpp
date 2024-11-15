#include <gtest/gtest.h>
#include <bcli/bcli.hpp>

using namespace bc;

TEST(Parser, simple)
{
  char* argv[] = {"cmd", "-p", "10", "pos1", "-t", "strvalue", "-f", "pos2", "pos3"};
  int argc = sizeof(argv)/sizeof(char*);
  {
    Parser cli("test", "test", "test", "test");

    cli.add_param("-p/--param", "help");
    cli.add_param("-t", "help");
    cli.add_param("-f", "help")->as_flag();

    BCLI_PARSE(cli, argc, argv)

    EXPECT_EQ(cli.get_positionals()[0], "pos1");
    EXPECT_EQ(cli.get_positionals()[1], "pos2");
    EXPECT_EQ(cli.get_positionals()[2], "pos3");

    int p = cli.getp("p")->as<int>();
    EXPECT_EQ(p, 10);

    EXPECT_EQ(cli.getp("t")->as<std::string>(), "strvalue");

    EXPECT_TRUE(cli.getp("f")->is_set());
  }
}

TEST(Parser, multiple)
{
  char* argv[] = {"cmd", "cmd1", "-p", "10", "pos1", "pos2", "pos3"};
  char* argv1[] = {"cmd", "cmd1", "-p", "10", "pos1", "pos2", "pos3", "pos4", "pos5", "pos6"};
  char* argv2[] = {"cmd", "cmd2", "-p", "42", "pos1", "pos2", "pos3"};
  int argc = sizeof(argv)/sizeof(char*);
  int argc1 = sizeof(argv1)/sizeof(char*);
  int argc2 = sizeof(argv2)/sizeof(char*);
  {
    Parser<1> cli("test", "test", "test", "test");
    cmd_t cmd1 = cli.add_command("cmd1", "cmd1 desc");
    cmd1->set_positional_bounds(1, 5, "Usage pos", "Help pos");

    cmd1->add_param("-p/--param", "help");
    cmd_t cmd2 = cli.add_command("cmd2", "cmd2 desc");

    EXPECT_THROW(cli.parse(argc1, argv1), ex::PositionalsError);
  }
  {
    Parser<1> cli("test", "test", "test", "test");
    
    cmd_t cmd1 = cli.add_command("cmd1", "cmd1 desc");
    cmd1->add_param("-p/--param", "help");

    cmd_t cmd2 = cli.add_command("cmd2", "cmd2 desc");
    cmd2->add_param("-p/--param", "help");

    cli.parse(argc2, argv2);
    EXPECT_EQ(cli.getp("p")->as<int>(), 42);

    cli.parse(argc, argv);
    EXPECT_EQ(cli.getp("p")->as<int>(), 10);

  }
}
