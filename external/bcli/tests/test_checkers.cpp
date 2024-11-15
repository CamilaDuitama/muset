#include <gtest/gtest.h>
#define _BCLI_TEST_
#include <bcli/bcli.hpp>

using namespace bc;

TEST(checkers, def)
{
  auto [b, s] = check::always_true("--test", "value");
  EXPECT_TRUE(b);
  EXPECT_EQ(s, "A true checker.");
}

TEST(checkers, is_file)
{
  {
    auto [res, msg] = check::is_file("--param", "./unknown.txt");
    EXPECT_FALSE(res);
  }
  {
    auto [res, msg] = check::is_file("--param", "./test_checkers.cpp");
    EXPECT_TRUE(res);
  }
}

TEST(checkers, is_dir)
{
  {
    auto [res, msg] = check::is_file("--param", "./unknown.dir");
    EXPECT_FALSE(res);
  }
  {
    auto [res, msg] = check::is_file("--param", "./");
    EXPECT_TRUE(res);
  }
}

TEST(checkers, is_number)
{
  {
    auto [res, msg] = check::is_number("--param", "123456789A");
    EXPECT_FALSE(res);
  }
  {
    auto [res, msg] = check::is_number("--param", "123456789");
    EXPECT_TRUE(res);
  }
}

TEST(checkers, is_dna)
{
  {
    auto [res, msg] = check::is_dna("--param", "ACGATTCGACGA");
    EXPECT_TRUE(res);
  }
  {
    auto [res, msg] = check::is_dna("--param", "ACGAXACGACGA");
    EXPECT_FALSE(res);
  }
}

TEST(checkers, is_rna)
{
  {
    auto [res, msg] = check::is_rna("--param", "ACGAUUUCGACGA");
    EXPECT_TRUE(res);
  }
  {
    auto [res, msg] = check::is_rna("--param", "ACGATTACGACGA");
    EXPECT_FALSE(res);
  }
}

TEST(checkers, ext)
{
  std::string fa = "./a/path/to/file.fasta";
  std::string fn = "./a/path/to/file.fna";
  std::string f = "./a/path/to/file.fa";
  std::string fql = "./a/path/to/file.fastq";
  std::string fq = "./a/path/to/file.fq";

  EXPECT_TRUE(std::get<0>(check::seems_fasta("--param", fa)));
  EXPECT_TRUE(std::get<0>(check::seems_fasta("--param", fn)));
  EXPECT_TRUE(std::get<0>(check::seems_fasta("--param", f)));
  
  EXPECT_TRUE(std::get<0>(check::seems_fastq("--param", fql)));
  EXPECT_TRUE(std::get<0>(check::seems_fastq("--param", fq)));

  EXPECT_FALSE(std::get<0>(check::seems_fasta("--param", fq)));
  EXPECT_FALSE(std::get<0>(check::seems_fastq("--param", fa)));
}

TEST(checkers, in)
{
  EXPECT_TRUE(std::get<0>(check::f::in("all|test|abc")("--param", "all")));
  EXPECT_FALSE(std::get<0>(check::f::in("alls|test|abc")("--param", "all")));
}

TEST(checkers, range)
{
  std::string r = "100";
  EXPECT_TRUE(std::get<0>(check::f::range(50, 100)("--param", r)));
  EXPECT_TRUE(std::get<0>(check::f::range(100, 101)("--param", r)));
  EXPECT_FALSE(std::get<0>(check::f::range(50, 99)("--param", r)));
  EXPECT_FALSE(std::get<0>(check::f::range(101, 150)("--param", r)));

  std::string d = "0.18";
  EXPECT_TRUE(std::get<0>(check::f::range(0.17, 0.20)("--param", d)));
  EXPECT_FALSE(std::get<0>(check::f::range(0.23, 0.99)("--param", d)));
}

TEST(checkers, lh)
{
  std::string r = "100";
  EXPECT_TRUE(std::get<0>(check::f::lower(101)("--param", r)));
  EXPECT_FALSE(std::get<0>(check::f::lower(99)("--param", r)));
  EXPECT_FALSE(std::get<0>(check::f::higher(101)("--param", r)));
  EXPECT_TRUE(std::get<0>(check::f::higher(99)("--param", r)));
}


TEST(checkers, check_magic)
{
  EXPECT_TRUE(std::get<0>(check::is_gz("--gz", "./data/test.txt.gz")));
  EXPECT_TRUE(std::get<0>(check::is_lz4_frame("--lz4", "./data/test.txt.lz4")));
  EXPECT_TRUE(std::get<0>(check::is_bz2("--bz2", "./data/test.txt.bz2")));
}