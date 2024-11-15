#include <gtest/gtest.h>
#define _BCLI_TEST_
#include <bcli/bcli.hpp>

using namespace bc;

TEST(param, make)
{
  param::param_t p = param::make("-p/--param", "make test");
  EXPECT_EQ(p->m_raw_name, "-p/--param");
  EXPECT_EQ(p->m_help, "make test");
  EXPECT_EQ(p->m_short, "p");
  EXPECT_EQ(p->m_long, "param");

  param::make("/p", "invalid param");
  EXPECT_THROW(ex::ExHandler::get().throw_last(), ex::BCliError);
  ex::ExHandler::get().clear();
}

TEST(param, param_setter)
{
  {
    std::string value;
    auto setter = [&value](const std::string& s) {value = s;};
    param::param_t p = param::make("-p/--param", "make test", setter);
    p->process("a random value");
    EXPECT_EQ(value, "a random value");
  }
  {
    int value;
    auto setter = [&value](const std::string& s) {value = std::stoi(s);};
    param::param_t p = param::make("-p/--param", "make test", setter);
    p->process("10");
    EXPECT_EQ(value, 10);
  }
  {
    int value;
    param::param_t p = param::make("-p/--param", "make test");
    p->setter(value);
    p->process("10");
    EXPECT_EQ(value, 10);

  }
}

TEST(param, param_checker)
{
  auto checker = check::is_number;
  {
    param::param_t p = param::make("-p/--param", "make test", nullptr, checker);
    p->process("10");
    EXPECT_TRUE(p->m_has_valid_value);
  }
  {
    param::param_t p = param::make("-p/--param", "make test", nullptr, checker);
    EXPECT_THROW(p->process("ZZ"), ex::CheckFailedError);
  }
  {
    param::param_t p = param::make("-p/--param", "make test");
    p->checker(check::is_number);
    EXPECT_THROW(p->process("ZZ"), ex::CheckFailedError);
  }
  {
    param::param_t p = param::make("-p/--param", "make test");
    p->checker(check::is_number);
    EXPECT_NO_THROW(p->process("10"));
    EXPECT_TRUE(p->m_has_valid_value);
  }
}

TEST(param, param_checker_setter)
{
  int value;
  auto checker = check::is_number;
  auto setter = param::get_setter(value);
  {
    param::param_t p = param::make("-p/-param", "make test");
    p->setter(setter)->checker(checker);
    EXPECT_NO_THROW(p->process("42"));
    EXPECT_EQ(value, 42);

  }
  {
    param::param_t p = param::make("-p/-param", "make test");
    p->setter(setter)->checker(checker);
    EXPECT_THROW(p->process("ZZ"), ex::CheckFailedError);
  }
}

TEST(param, param_group)
{
  param::pgroup_t pg = param::make_group("name", "desc");
  EXPECT_EQ(pg->m_name, "name");
  EXPECT_EQ(pg->m_desc, "desc");

  pg->add(param::make("-p/--param", "make test"));

  pg->add(param::make("-p/--other", "make test"));
  EXPECT_THROW(ex::ExHandler::get().check(), ex::BCliError);
  
  pg->add(param::make("-o/--param", "make test"));
  EXPECT_THROW(ex::ExHandler::get().check(), ex::BCliError);


  EXPECT_NO_THROW(pg->get("p"));
  EXPECT_NO_THROW(pg->get("param"));
}

TEST(param, command)
{
  param::cmd_t cmd = param::make_cmd("command", "description command");
  
  param::pgroup_t pg = param::make_group("main group", "main param");
  param::pgroup_t pg2 = param::make_group("advanced", "advanced params");

  cmd->add(pg);
  cmd->add(pg);
  EXPECT_THROW(ex::ExHandler::get().throw_last(), ex::BCliError);
  
  EXPECT_THROW(cmd->get("unknown group"), ex::UnknownGroupError);
  
  ex::ExHandler::get().clear();
  cmd->add(pg2);
  EXPECT_NO_THROW(ex::ExHandler::get().throw_last());
  
  cmd->get("advanced");
  EXPECT_NO_THROW(ex::ExHandler::get().throw_last());
}