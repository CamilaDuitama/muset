#include <vector>
#include <string>
#include <unordered_set>
#include <gtest/gtest.h>
#define _BCLI_TEST_
#include <bcli/bcli.hpp>

using vstring = std::vector<std::string>;
using namespace bc;

class class_1
{
public:
  class_1(bool i) : i(i) {}
  bool i;
};

class class_2
{
public:
  class_2(const std::string& v) : v(v) {}
  std::string v;
};

class class_3
{
public:
  class_3(const class_2& class2) : class2(class2) {}
  class_2 class2;
};

TEST(utils, lexical_cast)
{
  std::string s = "42";
  std::string s2 = "42.42";
  std::string s3 = "1e-10";
  int v = 42;
  double v2 = 42.42;
  double v3 = 1e-10;
  EXPECT_EQ(utils::lexical_cast<std::string>(s), s);
  EXPECT_EQ(utils::lexical_cast<int>(s), v);
  EXPECT_EQ(utils::lexical_cast<int>("42"), 42);
  EXPECT_EQ(utils::lexical_cast<double>(s2), v2);
  EXPECT_EQ(utils::lexical_cast<std::string>(v), s);
  EXPECT_EQ(utils::lexical_cast<std::string>(v2), s2);
  EXPECT_EQ(utils::lexical_cast<double>(s3), v3);
  EXPECT_EQ(utils::lexical_cast<std::string>(v3), s3);
  EXPECT_EQ(utils::lexical_cast<bool>(1), true);
  EXPECT_EQ(utils::lexical_cast<bool>(0), false);
  EXPECT_EQ(utils::lexical_cast<bool>(-1), true);
  EXPECT_EQ(utils::lexical_cast<bool>("1"), true);
  EXPECT_EQ(utils::lexical_cast<bool>("0"), false);
  EXPECT_EQ(utils::lexical_cast<std::string>(false), "0");
  EXPECT_EQ(utils::lexical_cast<std::string>(true), "1");
  EXPECT_EQ(utils::lexical_cast<class_1>(true).i, true);
  EXPECT_EQ(utils::lexical_cast<class_2>("TEST").v, "TEST");
  EXPECT_EQ(utils::lexical_cast<class_3>(class_2{"AAAA"}).class2.v, "AAAA");
  EXPECT_THROW(utils::lexical_cast<int>("."), ex::LexicalCastError);
}

TEST(utils, wrap)
{
  std::string s = "WRAP";
  std::string s2 = "[WRAP]";
  std::string s3 = "([WRAP])";
  EXPECT_EQ(utils::wrap(s, "[]"), s2);
  EXPECT_EQ(utils::wrap(s, "([])"), s3);

  EXPECT_EQ(utils::unwrap(s2, 1), s);
  EXPECT_EQ(utils::unwrap(s3, 2), s);
}

TEST(utils, join)
{
  vstring v = {"a", "b"};
  std::string del = "|";
  EXPECT_EQ(utils::join(v, "|"), "a|b");

  auto format = [](const std::string& s) -> std::string {
    return std::string("(") + s + std::string(")");
  };

  EXPECT_EQ(utils::join(v, "|", format), "(a)|(b)");
}

TEST(utils, contains)
{
  std::string s = "AAACGCTCTTTTACGACGA";
  EXPECT_TRUE(utils::contains(s, "TTTT"));
  EXPECT_FALSE(utils::contains(s, "AAAA"));
}
TEST(utils, start_end)
{
  std::string v = "ACGCATCACAGCA";
  EXPECT_TRUE(utils::startswith(v, "ACGCA"));
  EXPECT_FALSE(utils::startswith(v, "CCGCA"));
  EXPECT_TRUE(utils::endswith(v, "AGCA"));
  EXPECT_FALSE(utils::endswith(v, "CCGCA"));
}

TEST(utils, trim)
{
  std::string left = "  left";
  std::string right = "right   ";
  std::string both = "  both      ";

  EXPECT_EQ(utils::ltrim(left), "left");
  EXPECT_EQ(utils::rtrim(right), "right");
  EXPECT_EQ(utils::trim(both), "both");
}

TEST(utils, split)
{
  std::string s = "-t|--test";
  vstring v     = {"-t", "--test"};
  vstring vf    = {"(-t)", "(--test)"};
  vstring sp    = utils::split(s, '|');

  for (size_t i=0; i<v.size(); i++)
    EXPECT_EQ(v[i], sp[i]);

  auto format = [](const std::string& s) -> std::string {
    return std::string("(") + s + std::string(")");
  };

  vstring spf   = utils::split(s, '|', format);
  
  for (size_t i=0; i<vf.size(); i++)
    EXPECT_EQ(vf[i], spf[i]);

  std::string s2 = "test";
  vstring     v2 = {"te", "st"};
  auto [left, right] = utils::split2(s2);
  EXPECT_EQ(v2[0], left);
  EXPECT_EQ(v2[1], right);
}

TEST(utils, params)
{
  std::string lp = "--test";
  std::string sp = "-t";

  EXPECT_TRUE(utils::is_param(sp));
  EXPECT_TRUE(utils::is_param(lp));
  EXPECT_TRUE(utils::is_long_param(lp));
  EXPECT_TRUE(utils::is_short_param(sp));
  EXPECT_FALSE(utils::is_long_param(sp));
  EXPECT_FALSE(utils::is_short_param(lp));
}

TEST(utils, format)
{
  std::string f = utils::format_error("--param", "value", "an error msg");
  EXPECT_EQ(f, "[--param value] ~ an error msg");
}

TEST(utils, valid)
{
  auto validator = [](const char& c) -> bool {
    switch (c)
    {
      case 'A': return true;
      case 'C': return true;
      case 'T': return true;
      case 'G': return true;
      default: return false;
    }
  };

  std::string valid_dna = "ACGACTACTATCATA";
  std::string non_valid_dna = "AGAGACGATLAA";

  EXPECT_TRUE(utils::valid_string(valid_dna, validator));
  EXPECT_FALSE(utils::valid_string(non_valid_dna, validator));

}

TEST(utils, trim_params)
{
  std::string lp = "--test";
  std::string sp = "-t";

  EXPECT_EQ(utils::trim_param(sp), "t");
  EXPECT_EQ(utils::trim_param(lp), "test");
}