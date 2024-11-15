#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <filesystem>
#include <cstdio>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <unistd.h>
#include <fcntl.h>
#include <sys/sendfile.h>

#include <ext/stdio_filebuf.h>

#ifdef _CFRCAT_APP_
#include "options.hpp"
#endif

namespace fs = std::filesystem;

namespace cfr
{

using jumper_t = std::function<uint64_t(int)>;

template<typename... Path>
void fexist(Path... paths);

template<typename P, typename... Path>
void fexist(P path, Path... paths)
{
  if (!fs::exists(path))
    throw std::runtime_error(path + std::string(" doesn't exist!"));
  fexist(paths...);
}

template<>
void fexist() {return;}

int fileno_hack(const std::fstream& file_stream)
{
  return static_cast<__gnu_cxx::stdio_filebuf<char>*const>(file_stream.rdbuf())->fd();
}

FILE* cfile_hack(const std::fstream& file_stream)
{
  return static_cast<__gnu_cxx::stdio_filebuf<char>*const>(file_stream.rdbuf())->file();
}

int safe_open(const std::string& path, int flag, mode_t mode = 0, bool bypass = false)
{
  if (!bypass) fexist(path);

  int fd;
  if (mode)
    fd = open(path.c_str(), flag, mode);
  else
    fd = open(path.c_str(), flag);
  if (fd < 0)
    throw std::runtime_error("");
  return fd;
}


std::vector<int> open_files(const std::vector<std::string>& paths, int flag)
{
  std::vector<int> fds;
  for (auto& path: paths)
    fds.push_back(safe_open(path, flag));
  return fds;
}

void close_files(const std::vector<int>& fds)
{
  for (auto& fd: fds)
    close(fd);
}

uint64_t sizeof_fd(int fd)
{
  uint64_t size = lseek64(fd, 0, SEEK_END);
  lseek64(fd, 0, SEEK_SET);
  return size;
}

uint64_t walk_n(int fd, int n)
{
  return lseek64(fd, n, SEEK_SET);
}

uint64_t walk_until(int fd, char value, int n)
{
  if (!n) return 0;
  char buf[4096];
  uint64_t offset = 0, count = 0;

  for (int rd; rd=read(fd, buf, sizeof(buf));)
  {
    for (int i=0; i<rd; i++)
    {
      if (buf[i] != value)
        offset++;
      else
        count++;

      if (count == n) goto found;
    }
  }
  found:
    return lseek64(fd, offset+1, SEEK_SET);
  return 0;
}

uint64_t walk_n_lines(int fd, int n)
{
  return walk_until(fd, '\n', n);
}

uint64_t concat(int in_fd, int out_fd)
{
  int64_t offset = lseek64(out_fd, 0, SEEK_END);
  uint64_t size = sizeof_fd(in_fd);
  return copy_file_range(in_fd, NULL, out_fd, &offset, size, 0);
}

uint64_t concat(int in_fd, int out_fd, jumper_t jumper)
{
  int64_t offset = lseek64(out_fd, 0, SEEK_END);
  uint64_t size = sizeof_fd(in_fd);
  jumper(in_fd);
  return copy_file_range(in_fd, NULL, out_fd, &offset, size, 0);
}

void concat(int out_fd, std::vector<int>& fds, jumper_t jumper)
{
  for (auto& fd: fds)
    concat(fd, out_fd, jumper);
}

void concat(int out_fd, std::vector<int>& fds, std::vector<jumper_t>& jumpers)
{
  for (size_t i=0; i<fds.size(); i++)
    concat(fds[i], out_fd, jumpers[i]);
}

void concat(int fd_out, std::vector<std::string>& ins, jumper_t jumper)
{
  std::vector<int> fds = open_files(ins, O_RDONLY);
  concat(fd_out, fds, jumper);
  close_files(fds);
}

void concat(const std::string& out, std::vector<std::string>& ins, jumper_t jumper)
{
  int fd_out = safe_open(out, O_CREAT | O_TRUNC | O_RDWR, S_IRWXU, true);
  std::vector<int> fds = open_files(ins, O_RDONLY);
  concat(fd_out, fds, jumper);
  close_files(fds);
}

void copy_file(const std::string& in, const std::string& out)
{
  fexist(in, out);
  int in_fd = safe_open(in, O_RDONLY);
  int out_fd = safe_open(out, O_CREAT | O_TRUNC | O_RDWR, 0x01B6, true);
  uint64_t size = sizeof_fd(in_fd);
  sendfile(out_fd, in_fd, NULL, size);
}

void copy_file(int in, int out)
{
  uint64_t size = sizeof_fd(in);
  sendfile(out, in, NULL, size);
}

#ifdef _CFRCAT_APP_
void process(cfrcat_params& params)
{
  jumper_t jumper;
  if (params.skip_bytes)
    jumper = std::bind(&walk_n, std::placeholders::_1, params.skip_bytes);
  else if (params.skip_lines)
    jumper = std::bind(&walk_n_lines, std::placeholders::_1, params.skip_lines);
  else if (params.skip_char && params.n)
    jumper = std::bind(&walk_until, std::placeholders::_1, params.skip_char, params.n);
  else jumper = [](int fd) -> uint64_t {return lseek64(fd, 0, SEEK_SET);};
  
  int fd_out;
  if (params.out == "@")
    fd_out = STDOUT_FILENO;
  else if (!fs::exists(params.out) || !params.append)
    fd_out = safe_open(params.out, O_CREAT | O_TRUNC | O_RDWR, 0x01B6, true);
  else
    fd_out = safe_open(params.out, O_RDWR, 0, true);
  
  if (!params.header.empty())
  {
    int fd_header = safe_open(params.header, O_RDONLY);
    concat(fd_header, fd_out);
    close(fd_header);
  }
  concat(fd_out, params.ins, jumper);

  close(fd_out);
}
#endif
}