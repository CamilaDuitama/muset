/** MIT License
 *
 * Copyright (c) 2021 Téo Lemane
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * @file bcli.hpp
 * @author Téo Lemane
 * @version 0.0.1
*/

#pragma once
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <functional>
#include <filesystem>
#include <charconv>
#include <iomanip>
#include <vector>
#include <array>
#include <unordered_map>

#include <cassert>

/**
 * @mainpage
 *
 * Public API is described in Modules section. \n \n
 *
 * Modules:
 *  - @link Configuration @endlink
 *  - @link Exceptions @endlink
 *  - @link Utilities @endlink
 *  - @link Checkers @endlink
 *  - @link Param @endlink
 *  - @link ParamGroup @endlink
 *  - @link Command @endlink
 *  - @link Parser @endlink
 *
 * Examples:
 *  - @link 01_bcli_basic.cpp @endlink
 *  - @link 02_bcli_setter.cpp @endlink
 *  - @link 03_bcli_checker.cpp @endlink
 *  - @link 04_bcli_deps_bans.cpp @endlink
 *  - @link 05_bcli_get.cpp @endlink
 *  - @link 06_bcli_multiple.cpp @endlink
 *  - @link 07_bcli_real.cpp @endlink
 *  - @link 08_bcli_help.cpp @endlink
 *  - @link 09_bcli_help_multiple.cpp @endlink
 */

/**
 * @example 01_bcli_basic.cpp
 * @brief An example of basic usage.
 */

/**
 * @example 02_bcli_setter.cpp
 * @brief An example of setter usage.
 */

/**
 * @example 03_bcli_checker.cpp
 * @brief An example of checker usage.
 */

/**
 * @example 04_bcli_deps_bans.cpp
 * @brief An example of bcli param dependencies.
 */

/**
 * @example 05_bcli_get.cpp
 * @brief An example of bcli get param.
 */

/**
 * @example 06_bcli_multiple.cpp
 * @brief An example of bcli usage in multiple commands mode.
 */

/**
 * @example 07_bcli_real.cpp
 * @brief An example of real use case.
 */

/**
 * @example 08_bcli_help.cpp
 * @brief An example of customized help.
 */

/**
 * @example 09_bcli_help_multiple.cpp
 * @brief An example of customized help in multiple commands mode.
 */

#define BCLI_VERSION "0.0.1"

#ifdef _BCLI_TEST_
#define PRIVATE public
#else
#define PRIVATE private
#endif

/**
 * @namespace bc
 * @brief bcli namespace
 *
 * Main namespace.
 */
namespace bc {

namespace fs = std::filesystem;

template<typename... Args>
void unused(Args&&...) {}

/**
 * @defgroup Configuration
 * @brief About bcli configuration
 */

template<int Mode>
class Parser;

/**
 * @ingroup Configuration
 * @namespace config
 * @brief bcli configuration namespace
 *
 * This namespace contains essentially a config singleton (bc::config::Config).
 */
namespace config {

class Command;
class Param;

/**
 * @ingroup Configuration
 * @brief A configuration singleton to manage bcli configuration
 *
 * Shortcut bc::conf
 *
 * Usage:
 * @code
 * bc::conf::get().help(true)
 *                .verbose(true)
 *                .debug(false)
 *                .version(true)
 *                .flag_symbol("[F]")
 *                .default_group("main")
 *                .default_meta("STR");
 * @endcode
 */
class Config
{

private:
  Config() {}

public:
  /**
   * @brief Get the config singleton
   *
   * @return Config&
   */
  static Config& get()
  {
    static Config m_singleton;
    return m_singleton;
  }

  Config& help(bool v) {m_help = v; return *this;}
  Config& verbose(bool v) {m_verbose = v; return *this;}
  Config& debug(bool v) {m_debug = v; return *this;}
  Config& version(bool v) {m_version = v; return *this;}
  Config& flag_symbol(const std::string& sf) {m_flag_symbol = sf; return *this;}
  Config& default_group(const std::string& dg) {m_default_grp = dg; return *this;}
  Config& default_meta(const std::string& meta) {m_default_meta = meta; return *this;}

  bool has_common() {return m_help || m_verbose || m_debug || m_version;}

public:
  bool m_help {true};
  bool m_verbose {true};
  bool m_debug {true};
  bool m_version {true};

  std::string m_default_grp {"global"};
  std::string m_flag_symbol {"⚑"};
  std::string m_default_meta {"?"};

  friend class Parser<0>;
  friend class Parser<1>;
  friend class Param;
  friend class Command;
};

} // end of namespace config


/**
 * @defgroup Exceptions
 * @brief About bcli exceptions.
 *
 * Exceptions are contained within two groups:
 * - ImplError: handled by bc::ex::ExHandler
 *   - bc::ex::InvalidParamError
 *   - bc::ex::UnknownParamError
 *   - bc::ex::UnknownGroupError
 *   - bc::ex::CmdModeError
 *   - bc::ex::AlreadyExistsError
 * - UsageError: handled by the main try/catch
 *   - bc::ex::FileNotFoundError
 *   - bc::ex::DirNotFoundError
 *   - bc::ex::CheckFailedError
 *   - bc::ex::MissingValueError
 *   - bc::ex::UnknownCmdError
 *   - bc::ex::RequiredParamError
 *   - bc::ex::IncompatibleError
 *   - bc::ex::BannedError
 *   - bc::ex::DependsError
 *   - bc::ex::PostionalsError
 */

/**
 * @namespace ex
 * @ingroup Exceptions
 * @brief bcli exception namespace
 *
 * This namespace contains essentially bcli exceptions
 */
namespace ex {

/**
 * @enum ExitCodes
 * @ingroup Exceptions
 * @brief bcli exit codes
 *
 */
enum ExitCodes
{
  Sucess = 0,       /*!< Success. */
  Failure = 1,      /*!< Failure. */
  ImplError = 2,    /*!< Implementation error, before parsing. */
  UsageError = 3,   /*!< Usage error, during parsing. */
  UnknownError = 4  /*!< Unknow error. */
};

/**
 * @class BCliError
 * @ingroup Exceptions
 * @brief bcli base error class
 *
 */
class BCliError : public std::exception
{
PRIVATE:
  std::string name {"BCliError"};
  std::string msg  {"Base error msg should never be printed"};
  ExitCodes   code {ExitCodes::Failure};

public:
  BCliError(const std::string& name, const std::string& msg, ExitCodes code)
    : name(name), msg(msg), code(code) {}
  std::string get_name() const { return name; }
  std::string get_msg () const { return msg; }
  ExitCodes   get_exit_code() const { return code; }

};

/**
 * @def ERROR_CTR
 * @ingroup Exceptions
 * @brief A shortcut for exception constructor
 *
 */
#define ERROR_CTR(name, code) name(const std::string& msg) : BCliError(#name, msg, code) {}

/**
 * @def ERROR_CLS
 * @ingroup Exceptions
 * @brief Define a new exception class that inherits from BCliError
 *
 * **Usage:**
 * @code
 * ERROR_CLS(MyExcept, ex::ExitCodes::Failure)
 * @endcode
 */
#define ERROR_CLS(name, code)    \
  class name : public BCliError  \
  {                              \
  public:                        \
    ERROR_CTR(name, code)        \
  };                             \

/**
 * @exception InvalidParamError
 * @ingroup Exceptions
 * @brief Thrown if param name is invalid.
 *
 */
ERROR_CLS(InvalidParamError, ExitCodes::ImplError)

/**
 * @exception UnknownParamError
 * @ingroup Exceptions
 * @brief Thrown if the used parameter is unknown.
 *
 */
ERROR_CLS(UnknownParamError, ExitCodes::ImplError)

/**
 * @exception UnknownGroupError
 * @ingroup Exceptions
 * @brief Thrown if the called group is unknown.
 *
 */
ERROR_CLS(UnknownGroupError, ExitCodes::ImplError)

/**
 * @exception CmdModeError
 * @ingroup Exceptions
 * @brief Unused.
 *
 */
ERROR_CLS(CmdModeError, ExitCodes::ImplError)

/**
 * @exception AlreadyExistsError
 * @ingroup Exceptions
 * @brief Thrown if a parameter already exists in the current command.
 *
 */
ERROR_CLS(AlreadyExistsError, ExitCodes::ImplError)

/**
 * @exception FileNotFoundError
 * @ingroup Exceptions
 * @brief Thrown if the parameter value is not a valid file.
 *
 */
ERROR_CLS(FileNotFoundError, ExitCodes::UsageError)

/**
 * @exception DirNotFoundError
 * @ingroup Exceptions
 * @brief Thrown if the parameter value is not a valid dir.
 *
 */
ERROR_CLS(DirNotFoundError, ExitCodes::UsageError)

/**
 * @exception CheckFailedError
 * @ingroup Exceptions
 * @brief Thrown is value checking failed.
 *
 */
ERROR_CLS(CheckFailedError, ExitCodes::UsageError)

/**
 * @exception MissingValueError
 * @ingroup Exceptions
 * @brief Thrown if a value-parameter is used without any value.
 *
 */
ERROR_CLS(MissingValueError, ExitCodes::UsageError)

/**
 * @exception UnknownCmdError
 * @ingroup Exceptions
 * @brief Thrown if a cmd is unknown.
 *
 */
ERROR_CLS(UnknownCmdError, ExitCodes::UsageError)

/**
 * @exception RequiredParamError
 * @ingroup Exceptions
 * @brief Thrown if a required parameter is not used.
 *
 */
ERROR_CLS(RequiredParamError, ExitCodes::UsageError)

/**
 * @exception IncompatibleError
 * @ingroup Exceptions
 * @brief Thrown if parameters are incompatibles.
 *
 */
ERROR_CLS(IncompatibleError, ExitCodes::UsageError)

/**
 * @exception BannedError
 * @ingroup Exceptions
 * @brief Thrown if values of several parameters are inconsistent.
 *
 */
ERROR_CLS(BannedError, ExitCodes::UsageError)

/**
 * @exception DependsError
 * @ingroup Exceptions
 * @brief Thrown if values of several parameters are mutually dependent.
 *
 */
ERROR_CLS(DependsError, ExitCodes::UsageError)


/**
 * @exception PositionalsError
 * @ingroup Exceptions
 * @brief Thrown if positionals are inconsistent.
 *
 */
ERROR_CLS(PositionalsError, ExitCodes::UsageError)

/**
 * @brief LexicalCastError
 * @ingroup Exceptions
 * @brief Thrown if lexical_cast<T> fails.
 */
ERROR_CLS(LexicalCastError, ExitCodes::Failure)

/**
 * @ingroup Exceptions
 * @brief bcli exception handler
 *
 * A singleton exception handler that aggregate all thrown exceptions about invalid impl
 * -> ex: invalid param name, already used param ect...
 *
 * ExHandler::check is called before parsing user args to ensure that cli is valid.
 * If not all exceptions are printed.
 *
 * Others exceptions are catched by the main try/catch, they concerns invalid usage and can be
 * throw during user args parsing.
 * -> ex: required param, incompatible params, value check failed etc...
 */
class ExHandler
{
public:
  static ExHandler& get()
  {
    static ExHandler m_singleton;
    return m_singleton;
  }

private:
  ExHandler() {}

public:
  void push(const BCliError& e)
  {
    m_exceptions.push_back(e);
  }

  void check()
  {
    if (!m_exceptions.empty())
    {
      for (auto& e : m_exceptions)
        std::cerr << "[" << e.get_name() << "]" << " -> " << e.get_msg() << std::endl;
      throw BCliError("", "", Failure);
    }
  }

  void clear()
  {
    m_exceptions = std::vector<BCliError>{};
  }

  void throw_last()
  {
    if (!m_exceptions.empty())
    {
      BCliError e = m_exceptions.back();
      m_exceptions.pop_back();
      throw e;
    }
  }

PRIVATE:
  std::vector<BCliError> m_exceptions;
};

} // end of namespace ex


/**
 * @defgroup Utilities
 * @brief About bcli utilities.
 */

/**
 * @ingroup Utilities
 * @namespace utils
 * @brief bcli utilities
 *
 */
namespace utils {

template<typename>
struct is_string : std::false_type {};

template<>
struct is_string<std::string> : std::true_type {};

template<>
struct is_string<std::string_view> : std::true_type {};

template<typename T>
constexpr auto is_string_v = is_string<T>::value;

template<typename T>
struct rwrapper {using type = T;};

template<>
struct rwrapper<std::string_view> {using type = std::string;};

template<typename RetType>
using rwrapper_t = typename rwrapper<RetType>::type;

enum class LexicalCast
{
  from_str,
  to_str,
  str,
  implicit,
  constructible
};

template<typename R, typename T, LexicalCast mode>
rwrapper_t<R> lexical_cast(const T& src)
{
  if constexpr(mode == LexicalCast::from_str)
  {
    std::stringstream ss;
    ss << src;
    if (R ret; ss >> ret)
      return ret;
    else
      throw ex::LexicalCastError(
        "Unable to cast \"" + std::string(src) + "\" to " + std::string(typeid(R).name()));
  }
  else if constexpr(mode == LexicalCast::to_str)
  {
    std::stringstream ss;
    ss << src;
    return ss.str();
  }
  else if constexpr(mode == LexicalCast::str)
  {
    return {std::begin(src), std::end(src)};
  }
  else if constexpr(mode == LexicalCast::constructible)
  {
    return R(src);
  }
  else if constexpr(mode == LexicalCast::implicit)
  {
    return src;
  }
}

template<typename R, typename T>
auto lexical_cast(T&& src)
{
  if constexpr(std::is_same_v<R, T> || std::is_convertible_v<R, T>)
    return lexical_cast<R, T, LexicalCast::implicit>(std::forward<T>(src));
  else if constexpr(std::is_constructible_v<R, T> && !std::is_same_v<R, bool>)
    return lexical_cast<R, T, LexicalCast::constructible>(std::forward<T>(src));
  else
    if constexpr(is_string_v<R>)
      if constexpr(is_string_v<T>)
        return lexical_cast<R, T, LexicalCast::str>(std::forward<T>(src));
      else
        return lexical_cast<R, T, LexicalCast::to_str>(std::forward<T>(src));
    else
      return lexical_cast<R, T, LexicalCast::from_str>(std::forward<T>(src));
}

/**
 * @ingroup Utilities
 * @brief join
 * @tparam T
 * @tparam Formatter=std::function<std::string(const std::string&)>
 * @param v an iterable string container
 * @param delim a string delimiter
 * @param format a formatter callback
 * @return std::string
 */
template<typename T, typename Formatter = std::function<std::string(const std::string&)>>
inline std::string join(const T& v, const std::string& delim, Formatter format)
{
  std::stringstream ss;
  for (auto it=std::begin(v); it!=std::end(v);)
    it != std::end(v)-1 ? ss << format(*it++) << delim : ss << format(*it++);
  return ss.str();
}

/**
 * @ingroup Utilities
 * @brief join
 *
 * @tparam T
 * @param v an iterable string container
 * @param delim a string delimiter
 * @return std::string
 */
template<typename T>
inline std::string join(const T& v, const std::string& delim)
{
  return join(v, delim, [](const std::string& s) -> std::string {return s;});
}

/**
 * @ingroup Utilities
 * @brief startswith
 * @param s
 * @param p
 * @return true if s starts with p
 * @return false if s doesn't start with p
 */
inline bool startswith(const std::string& s, const std::string& p)
{
  if (p.size() > s.size()) return false;
  return std::equal(p.begin(), p.end(), s.begin());
}

/**
 * @ingroup Utilities
 * @brief endswith
 * @param s
 * @param p
 * @return true if s ends with p
 * @return false if s doesn't end with p
 */
inline bool endswith(const std::string& s, const std::string& p)
{
  if (p.size()> s.size()) return false;
  return std::equal(p.rbegin(), p.rend(), s.rbegin());
}

/**
 * @ingroup Utilities
 * @brief contains
 *
 * @param s
 * @param p
 * @return true if s contains p
 * @return false if p not in s
 */
inline bool contains(const std::string& s, const std::string& p)
{
  if (s.find(p) != std::string::npos)
    return true;
  return false;
}

const std::string WHITECHAR = " \n\r\t\f\v";

/**
 * @ingroup Utilities
 * @brief rtrim
 *
 * Right trimming.
 * @code
 * trim("-test-|", "|-") -> "-test"
 * @endcode
 *
 * @param s a string to trim
 * @param v a string of banned char
 * @return std::string
 */
inline std::string rtrim(const std::string& s, const std::string& v = WHITECHAR)
{
  size_t end = s.find_last_not_of(v);
  if (end == std::string::npos)
    return "";
  return s.substr(0, end+1);
}

/**
 * @ingroup Utilities
 * @brief ltrim
 *
 * Left trimming.
 * @code
 * trim("-test-|", "|-") -> "test-|"
 * @endcode
 *
 * @param s a string to trim
 * @param v a string of banned char
 * @return std::string
 */
inline std::string ltrim(const std::string& s, const std::string& v = WHITECHAR)
{
  size_t beg = s.find_first_not_of(v);
  if (beg == std::string::npos)
    return "";
  return s.substr(beg);
}

/**
 * @ingroup Utilities
 * @brief trim
 *
 * Both sides trimming.
 * @code
 * trim("-test-|", "|-") -> "test"
 * @endcode
 *
 * @param s a string to trim
 * @param v a string of banned char
 * @return std::string
 */
inline std::string trim(const std::string& s, const std::string& v = WHITECHAR)
{
  return ltrim(rtrim(s, v), v);
}

/**
 * @ingroup Utilities
 *from_str_to_Tplit string in half.
 *
 * @code
 * split2("test") -> {"te", "st"}
 * @endcode
 *
 * @param s a string to split
 * @return std::tuple<std::string, std::string>
 */
inline std::tuple<std::string, std::string> split2(const std::string& s)
{
  assert(s.size() % 2 == 0);
  size_t mid = s.size() / 2;
  return std::make_tuple(s.substr(0, mid), s.substr(mid));
}

/**
 * @ingroup Utilities
 * @brief wrap
 *
 * @code
 * wrap("test", "[]") -> "[test]"
 * @endcode
 * @param s a string to wrap
 * @param bounds wrappers, see code example
 * @return std::string
 */
inline std::string wrap(const std::string& s, const std::string& bounds)
{
  std::stringstream ss;
  auto [l, r] = split2(bounds);
  ss << l << s << r;
  return ss.str();
}

/**
 * @ingroup Utilities
 * @brief unwrap
 *
 * @code
 * unwrap("([test])", 1) -> "[test]"
 * unwrap("([test])", 2) -> "test"
 * @endcode
 * @param s a string to unwrap
 * @param size size, see code example
 * @return std::string
 */
inline std::string unwrap(const std::string& s, size_t size)
{
  return s.substr(size, s.size()-2*size);
}

/**
 * @ingroup Utilities
 * @brief split
 *
 * @tparam std::function<std::string(const std::string&)>
 * @param s a string to split
 * @param delim a char delimiter
 * @param format a string formatter
 * @return std::vector<std::string>
 */
template<typename Formatter = std::function<std::string(const std::string&)>>
std::vector<std::string> split(const std::string& s, char delim, Formatter format)
{
  std::istringstream iss(s);
  std::vector<std::string> ret;
  for (std::string tmp; std::getline(iss, tmp, delim);)
    ret.push_back(format(tmp));
  return ret;
}

/**
 * @ingroup Utilities
 * @brief split
 *
 * @param s a string to split
 * @param delim a char delimiter
 * @return std::vector<std::string>
 */
inline std::vector<std::string> split(const std::string& s, char delim)
{
  return split(s, delim, [](const std::string& s) -> std::string {return s;});
}

/**
 * @ingroup Utilities
 * @brief is_long_param
 *
 * @code
 * is_long_param("--param") -> true
 * @endcode
 * @param p a parameter
 * @return true if p starts with "--"
 * @return false
 */
inline bool is_long_param(const std::string& p)
{
  return p.rfind("--", 0) == 0;
}

/**
 * @ingroup Utilities
 * @brief is_short_param
 *
 * @code
 * is_short_param("-p") -> true
 * @endcode
 * @param p a parameter
 * @return true if p starts with "-"
 * @return false
 */
inline bool is_short_param(const std::string& p)
{
  if (is_long_param(p))
    return false;
  return p.rfind("-", 0) == 0;
}

/**
 * @ingroup Utilities
 * @brief is_param
 *
 * @param p a parameter
 * @return true if p starts with "-" or "--"
 * @return false
 */
inline bool is_param(const std::string& p)
{
  return is_short_param(p) || is_long_param(p);
}

/**
 * @ingroup Utilities
 * @brief trim_param
 *
 * @code
 * trim_param("--param") -> "param"
 * @endcode
 * @param p a parameter to trim
 * @return std::string
 */
inline std::string trim_param(const std::string& p)
{
  return ltrim(p, "-");
}

/**
 * @ingroup Utilities
 * @brief format_error
 *
 * @param p
 * @param v
 * @param m
 * @return std::string
 */
inline std::string format_error(const std::string& p, const std::string& v, const std::string& m)
{
  std::stringstream ss;
  ss << "[" << p << " " << v << "]" << " ~ " << m;
  return ss.str();
}

/**
 * @ingroup Utilities
 * @brief  format_depend_errors
 *
 * @tparam T
 * @param p
 * @param d
 * @param msg
 * @return std::string
 */
template<typename T>
inline std::string format_depend_errors(T p, T d, const std::string& msg)
{
  std::stringstream ss;
  ss << p->raw() << " " << p->value() << " depends on " << d->raw();
  ss << " ~ " << "Check returns -> " << msg;
  return ss.str();
}

/**
 * @ingroup Utilities
 * @brief format_banned_errors
 *
 * @tparam T
 * @param p
 * @param d
 * @param msg
 * @return std::string
 */
template<typename T>
inline std::string format_banned_errors(T p, T d, const std::string& msg)
{
  std::stringstream ss;
  ss << p->raw() << " " << p->value() << " is banned with " << d->raw() << " " << d->value();
  ss << " ~ " << "Check returns -> " << msg;
  return ss.str();
}

class param_t;
template<>
inline std::string format_banned_errors(param_t, param_t, const std::string& msg);

/**
 * @ingroup Utilities
 * @brief valid_string
 *
 * @code
 * auto validator = [](const char& c) -> bool {
 *   if (c == 'A' || c == 'C' || c == 'T' || c == 'G')
 *     return true;
 *   return false;
 * };
 * valid_string("ACGTTTACGA", validator);
 * @endcode
 * @tparam std::function<bool(const char&)>
 * @param s
 * @param valid
 * @return true if s is valid
 * @return false if s is invalid
 */
template<typename Validator = std::function<bool(const char&)>>
bool valid_string(const std::string& s, Validator valid)
{
  for (const char& c : s)
    if (!valid(c))
      return false;
  return true;
}

/**
 * @ingroup Utilities
 * @brief sp
 *
 * @param n
 * @return std::string
 */
inline std::string sp(int n)
{
  return std::string(n, ' ');
}

/**
 * @ingroup Utilities
 * @brief exit_bcli
 *
 * @param e
 */
inline void exit_bcli(const ex::BCliError& e)
{
  if (!e.get_name().empty())
    std::cerr << wrap(e.get_name(), "[]") << " -> " << e.get_msg() << std::endl;
}

} // end of namespace utils

/**
 * @defgroup Checkers
 * @brief About bcli checkers.
 *
 * Checkers are functions used to parameter value validation.
 * A checker take two parameters:
 *  - an id, usually a parameter (ex: "--param")
 *  - a value to check
 *
 * And returns a std::tuple<bool, std::string>
 *  - 0: true if check success, false otherwise
 *  - 1: An error message
 *
 * Two namespaces:
 *  - check -> contains checkers
 *  - check::f -> contains factories
 *
 */

/**
 * @namespace check
 * @ingroup Checkers
 * @brief bcli check namespace
 *
 */
namespace check {

/**
 * @typedef checker_ret_t
 * @ingroup Checkers
 * @brief checkers return type
 *
 */
using checker_ret_t = std::tuple<bool, std::string>;

/**
 * @typedef checker_fn_t
 * @ingroup Checkers
 * @brief checkers function signature
 *
 */
using checker_fn_t = std::function<checker_ret_t(const std::string&, const std::string&)>;


/**
 * @ingroup Checkers
 * @brief throw ex::CheckFailedError if checker returns false
 *
 * ex: throw_if_false(check::is_number("--integer", "10"));
 */
inline void throw_if_false(const checker_ret_t& rc)
{
  if (!std::get<0>(rc))
    throw ex::CheckFailedError(std::get<1>(rc));
}

// Shortcut to define a new checker
/**
 * @ingroup Checkers
 * @def DEFINE_CHECKER
 * @brief define a new checker
 */
#define DEFINE_CHECKER(name, s1, s2)                                        \
  inline checker_ret_t name(const std::string& s1, const std::string& s2)   \


DEFINE_CHECKER(always_true, p, v)
{
  unused(p);
  unused(v);
  return std::make_tuple(true, "A true checker.");
}

/**
 * @ingroup Checkers
 * @brief is_file checker
 *
 * Return true if v corresponds to a file.
 *
 * ex: check::is_file("--param", "/path/to/file.txt");
 *
 * @param p parameter as string (ex: "--param")
 * @param v parameter value
 *
 * @return std::tuple<bool, std::string>
 */
DEFINE_CHECKER(is_file, p, v)
{
  bool exists = fs::exists(v);
  return std::make_tuple(exists, utils::format_error(p, v, "File doesn't exist!"));
}

/**
 * @ingroup Checkers
 * @brief is_dir checker
 *
 * Return true if v corresponds to a directory.
 *
 * ex: check::is_dir("--param", "/path/to/dir");
 *
 * @param p parameter as string (ex: "--param")
 * @param v parameter value
 *
 * @return std::tuple<bool, std::string>
 */
DEFINE_CHECKER(is_dir, p, v)
{
  bool exists = fs::is_directory(v);
  return std::make_tuple(exists, utils::format_error(p, v, "Directory doesn't exist!"));
}

/**
 * @ingroup Checkers
 * @brief is_number checker
 *
 * Return true if v corresponds to a number.
 *
 * ex: check::is_number("--param", "10");
 *
 * @param p parameter as string (ex: "--file")
 * @param v parameter value
 *
 * @return std::tuple<bool, std::string>
 */
DEFINE_CHECKER(is_number, p, v)
{
  bool is = true;
  bool neg = false;
  for (const char& c : v)
  {
    if (!neg && c == '-')
    {
      neg = true;
      continue;
    }
    if (!std::isdigit(c) && c != '.')
    {
      is = false;
      break;
    }
  }
  if (v.empty()) is = false;
  return std::make_tuple(is, utils::format_error(p, v, "Not a number!"));
}

/**
 * @ingroup Checkers
 * @brief is_dna checker
 *
 * Return true if v corresponds to a dna string.
 *
 * ex: check::is_dna("--param", "ACGACGA");
 *
 * @param p parameter as string (ex: "--param")
 * @param v parameter value
 *
 * @return std::tuple<bool, std::string>
 */
DEFINE_CHECKER(is_dna, p, v)
{
  auto valid_nt = [](const char& c) -> bool {
    switch (c)
    {
      case 'A': return true;
      case 'C': return true;
      case 'T': return true;
      case 'G': return true;
      default: return false;
    }
  };
  return std::make_tuple(utils::valid_string(v, valid_nt),
                         utils::format_error(p, v, "Not a valid dna string."));
}

/**
 * @ingroup Checkers
 * @brief is_number checker
 *
 * Return true if v corresponds to a rna string.
 *
 * ex: check::is_rna("--param", "10");
 *
 * @param p parameter as string (ex: "--param")
 * @param v parameter value
 *
 * @return std::tuple<bool, std::string>
 */
DEFINE_CHECKER(is_rna, p, v)
{
  auto valid_nt = [](const char& c) -> bool {
    switch (c)
    {
      case 'A': return true;
      case 'C': return true;
      case 'U': return true;
      case 'G': return true;
      default: return false;
    }
  };
  return std::make_tuple(utils::valid_string(v, valid_nt),
                         utils::format_error(p, v, "Not a valid rna string."));
}

/**
 * @namespace f
 * @ingroup Checkers
 * @brief bcli checker factories namespace
 *
 * Checker factories are function that return check::checker_fn_t.
 * This allows configurable checkers.
 *
 */
namespace f {

/**
 * @ingroup Checkers
 * @brief extension checker factory
 *
 * @code
 * auto seems_fasta = check::f::ext("fa|fna");
 * throw_if_false(seems_fasta("--param", "/path/to/file.fasta"));
 * @endcode
 * @param ext A set of extensions as string, sep by '|'.
 * @return checker_fn_t
 */
inline checker_fn_t ext(const std::string& ext)
{
  return [ext](const std::string& p, const std::string&v) -> checker_ret_t {
    std::vector<std::string> exts = utils::split(ext, '|');
    auto error_msg = std::bind(&utils::format_error, p, v, std::placeholders::_1);
    if (!fs::path(v).has_extension())
      return std::make_tuple(false, error_msg("No extension."));
    else
    {
      if (std::any_of(exts.begin(), exts.end(), std::bind(&utils::endswith, v, std::placeholders::_1)))
      {
        return std::make_tuple<bool, std::string>(true, {});
      }
      else
      {
        std::string extp = fs::path(v).extension();
        return std::make_tuple(false, error_msg(extp + "!=" + ext));
      }
    }
  };
}

/**
 * @ingroup Checkers
 * @brief range checker factory
 * @code
 * auto in_10_100 = check::range(10, 100);
 * throw_if_false(in_10_100("--param", "42"));
 * @endcode
 * @tparam T An arithmetic type
 * @tparam std::enable_if_t<std::is_arithmetic_v<T>, void>
 * @param start lower bound
 * @param end upper bound
 * @return checker_fn_t
 */
template<typename T,
         typename = typename std::enable_if_t<std::is_arithmetic_v<T>, void>>
inline checker_fn_t range(T start, T end)
{
  return [start, end](const std::string& p, const std::string& v) -> checker_ret_t {
    throw_if_false(is_number(p, v));
    T value;
    if constexpr(std::is_floating_point_v<T>)
      value = std::stod(v);
    else if (std::is_integral_v<T>)
      value = std::stoi(v);
    bool in_range = (value >= start) && (value <= end);
    return std::make_tuple(in_range, utils::format_error(
      p, v, "Not in range [" + std::to_string(start) + "," + std::to_string(end) +"]."));
  };
}

/**
 * @ingroup Checkers
 * @brief in checker factory
 * @code
 * auto is_valid_cmd = in("cmd2|cmd3");
 * throw_if_false(is_valid_cmd("--param", "cmd3"));
 * @endcode
 * @param s A set of strings as one string, sep by '|'.
 * @return checker_fn_t
 */
inline checker_fn_t in(const std::string& s)
{
  return [s](const std::string& p, const std::string& v) -> checker_ret_t {
    std::vector<std::string> vs = utils::split(s, '|');
    bool is_in = false;
    if (std::any_of(vs.begin(), vs.end(), [v](const std::string& vv) {
      return v == vv;
    })) {is_in = true;}

    return std::make_tuple(is_in, utils::format_error(
      p, v, "Not in " + utils::wrap(s, "[]")
    ));
  };
}

/**
 * @ingroup Checkers
 * @brief lower checker factory
 * @code
 * auto lower10 = lower(10);
 * throw_if_false(lower10("--param", "5"));
 * @endcode
 * @param n
 * @return checker_fn_t
 */
inline checker_fn_t lower(int n)
{
  return [n](const std::string& p, const std::string& v) -> checker_ret_t {
    throw_if_false(is_number(p, v));
    int value = std::stoi(v);
    return std::make_tuple(value < n, utils::format_error(p, v, v + ">=" + std::to_string(n)));
  };
}

/**
 * @ingroup Checkers
 * @brief higher checker factory
 * @code
 * auto higher10 = higher(10);
 * throw_if_false(lower10("--param", "15"));
 * @endcode
 * @param n
 * @return checker_fn_t
 */
inline checker_fn_t higher(int n)
{
  return [n](const std::string& p, const std::string& v) -> checker_ret_t {
    throw_if_false(is_number(p, v));
    int value = std::stoi(v);
    return std::make_tuple(value > n, utils::format_error(p, v, v + "<=" + std::to_string(n)));
  };
}

/**
 * @ingroup Checkers
 * @brief check_magic checker factory
 * @code
 * auto is_gz = f::check_magic<2>("gz", {0x1F, 0x8B});
 * throw_if_false(is_gz("--param", "/path/to/file.gz"));
 * @endcode
 * @tparam SIZE number of magic bytes, std::array<uint8_t, SIZE>
 * @param name a name
 * @param flag an array of bytes
 * @return checker_fn_t
 */
template<size_t SIZE>
inline checker_fn_t check_magic(const std::string& name, std::array<uint8_t, SIZE> flag)
{
  return [name, flag](const std::string& p, const std::string& v) -> checker_ret_t {
    throw_if_false(is_file(p, v));

    bool is_valid = false;

    std::array<uint8_t, SIZE> buffer;
    std::ifstream inf(v, std::ios::binary | std::ios::in);
    if (inf.good())
      inf.read(reinterpret_cast<char*>(buffer.data()), buffer.size());

    if (buffer == flag) is_valid = true;

    return std::make_tuple(is_valid,
                           utils::format_error(p, v, "Not a " + name + " file."));
  };
}

} // end of namespace f (checker factories)

/**
 * @ingroup Checkers
 * @brief fastx ext checker
 *
 */
inline checker_fn_t seems_fastx = f::ext("fa|fna|fasta|fastq|fq");

/**
 * @ingroup Checkers
 * @brief fasta ext checker
 *
 */
inline checker_fn_t seems_fasta = f::ext("fa|fna|fasta");

/**
 * @ingroup Checkers
 * @brief fastq ext checker
 *
 */
inline checker_fn_t seems_fastq = f::ext("fastq|fq");

/**
 * @ingroup Checkers
 * @brief sam ext checker
 *
 */
inline checker_fn_t seems_sam = f::ext("sam");

/**
 * @ingroup Checkers
 * @brief bam ext checker
 *
 */
inline checker_fn_t seems_bam = f::ext("bam");

/**
 * @ingroup Checkers
 * @brief cram ext checker
 *
 */
inline checker_fn_t seems_cram = f::ext("cram");

/**
 * @ingroup Checkers
 * @brief comp ext checker
 *
 */
inline checker_fn_t seems_comp = f::ext("gz|bz2|lz4");

/**
 * @ingroup Checkers
 * @brief gz ext checker
 *
 */
inline checker_fn_t seems_gz = f::ext("gz");

/**
 * @ingroup Checkers
 * @brief fastx
 *
 */
inline checker_fn_t seems_lz4 = f::ext("lz4");

/**
 * @ingroup Checkers
 * @brief gz checker
 *
 */
inline auto is_gz = f::check_magic<2>("gz", {0x1F, 0x8B});

/**
 * @ingroup Checkers
 * @brief lz4_frame checker
 *
 */
inline auto is_lz4_frame = f::check_magic<4>("lz4frame", {0x04, 0x22, 0x4D, 0x18});

/**
 * @ingroup Checkers
 * @brief bz2 checker
 *
 */
inline auto is_bz2 = f::check_magic<3>("bz2", {0x42, 0x5A, 0x68});

/**
 * @ingroup Checkers
 * @brief bam checker
 *
 */
inline auto is_bam = f::check_magic<4>("bam", {0x1F, 0x8B, 0x08, 0x04});

/**
 * @ingroup Checkers
 * @brief cram checker
 *
 */
inline auto is_cram = f::check_magic<4>("cram", {0x43, 0x52, 0x41, 0x4d});

} // end of namespace checker


/**
 * @defgroup Param
 * @brief About bcli parameters
 */

const std::string FLAG_VALUE = "BCLI_FSET";

/**
 * @ingroup Param
 * @enum Action
 * @brief parameter actions
 *
 */
enum Action
{
  Nothing,     /*!< Nothing  */
  ShowHelp,    /*!< Trigger help */
  ShowVersion  /*!< Trigger version */
};

/**
 * @namespace param
 * @brief bcli param namespace
 *
 *
 */
namespace param {

class Param;

/**
 * @ingroup Param
 * @typedef param_t
 * @brief param_t
 *
 * Param constuctor is private, Param is always used as std::shared_ptr<Param>
 */
using param_t = std::shared_ptr<Param>;

/**
 * @ingroup Param
 * @typedef checker_fn_t
 * @brief checker function signature
 *
 *
 */
using checker_fn_t = check::checker_fn_t;

/**
 * @ingroup Param
 * @typedef setter_fn_t
 * @brief setter function signature
 *
 */
using setter_fn_t = std::function<void(const std::string&)>;

/**
 * @ingroup Param
 * @typedef help_fn_t
 * @brief help function signature
 *
 */
using help_fn_t = std::function<std::string()>;

/**
 * @ingroup Param
 * @typedef callback_fn_t
 * @brief callback function signature
 *
 */
using callback_fn_t = std::function<void()>;

using conf = config::Config;

/**
 * @brief get a setter from variable reference
 *
 * Returns a setter to set T var
 *
 * @tparam T
 * @param var variable to set
 * @return setter_fn_t
 */
template<typename T>
setter_fn_t get_setter(T& var)
{
  return [&var] (const std::string& in = {}) {
    if constexpr(std::is_same_v<T, bool>)
    {
      if (in == FLAG_VALUE)
        var = true;
      else
        var = false;
    }
    else
    {
      var = utils::lexical_cast<T>(in);
    }
  };
}

//enum CheckerMode
//{
//  AND,
//  OR,
//  XOR
//};

class ParamGroup;

/**
 * @ingroup Param
 * @brief Param
 *
 * Parameter class
 */
class Param: public std::enable_shared_from_this<Param>
{
  using s_ptr = std::shared_ptr<Param>;
  friend param_t make(const std::string& name,
                      const std::string& help,
                      setter_fn_t setter,
                      checker_fn_t checker);

  friend class ParamTest;
  friend class ParamGroup;
  friend class Command;
  friend class Parser<0>;
  friend class Parser<1>;

public:
  /**
   * @brief get idx
   *
   * @return std::tuple<std::string, std::string>
   */
  std::tuple<std::string, std::string> idx()
  {
    return std::make_tuple(m_short, m_long);
  }

  /**
   * @brief get raw name
   *
   * @return std::string
   */
  std::string raw()
  {
    return m_raw_name;
  }

  /**
   * @brief set default value
   *
   * @param default_value
   * @return param_t
   */
  param_t def(const std::string& default_value)
  {
    m_default = default_value;
    m_str_value = m_default;
    m_has_default = true;
    return shared_from_this();
  }

  /**
   * @brief set checker
   *
   * @param checker_callback
   * @return param_t
   */
  param_t checker(checker_fn_t checker_callback)
  {
    if (m_is_flag)
    {
      ex::ExHandler::get().push(
        ex::IncompatibleError(utils::wrap(m_raw_name, "[]") + "~ A flag cannot have a checker."));
    }
    c_checkers.push_back(checker_callback);
    return shared_from_this();
  }

  //param_t checker_mode(CheckerMode mode)
  //{
  //  m_check_mode = mode;
  //}

  /**
   * @brief set setter (auto from variable reference)
   *
   * @tparam T
   * @param var
   * @return param_t
   */
  template<typename T>
  param_t setter(T& var)
  {
    c_setter = get_setter<T>(var);
    return shared_from_this();
  }

  /**
   * @brief set setter
   *
   * @param setter_callback
   * @return param_t
   */
  param_t setter_c(setter_fn_t setter_callback)
  {
    c_setter = setter_callback;
    return shared_from_this();
  }

  /**
   * @brief set callback
   *
   * Callback is called if parameter is used by user.
   *
   * @param callback
   * @return param_t
   */
  param_t callback(callback_fn_t callback)
  {
    c_callback = callback;
    m_callback_trigger = true;
    return shared_from_this();
  }

  /**
   * @brief set dependencies
   *
   * see @link 04_bcli_deps_bans.cpp @endlink
   *
   * @param checker
   * @param p
   * @param pchecker
   * @return param_t
   */
  param_t depends_on(check::checker_fn_t checker, param_t p, check::checker_fn_t pchecker = nullptr)
  {
    m_depends_on.push_back({checker, p, pchecker});
    return shared_from_this();
  }

  /**
   * @brief
   *
   * see @link 04_bcli_deps_bans.cpp @endlink
   *
   * @param checker
   * @param p
   * @param pchecker
   * @return param_t
   */
  param_t banned(check::checker_fn_t checker, param_t p, check::checker_fn_t pchecker = nullptr)
  {
    m_banned.push_back({checker, p, pchecker});
    return shared_from_this();
  }

  /**
   * @brief use param as flag (without value)
   *
   * @return param_t
   */
  param_t as_flag()
  {
    if (c_checker)
    {
      ex::ExHandler::get().push(
        ex::IncompatibleError(
          utils::wrap(m_raw_name, "[]") + "~ A param with checker cannot be a flag."));
    }
    m_is_flag = true;
    m_has_default = true;
    return shared_from_this();
  }

  /**
   * @brief set metavar
   *
   * Metavar shown in usage: \n
   * meta("INT") -> -i/--integer <INT> \n
   * meta("FILE") -> --file <FILE> \n
   *
   * @param meta
   * @return param_t
   */
  param_t meta(const std::string& meta)
  {
    m_meta = meta;
    return shared_from_this();
  }

  /**
   * @brief set action
   *
   * @see bc::Action
   *
   * @param action
   * @return param_t
   */
  param_t action(Action action)
  {
    m_action = action;
    return shared_from_this();
  }

  param_t hide()
  {
    m_hidden = true;
    return shared_from_this();
  }

  /**
   * @brief get str value
   *
   * @return const std::string&
   */
  const std::string& value()
  {
    return m_str_value;
  }

  /**
   * @brief get default value
   *
   * @return const std::string&
   */
  const std::string& get_def()
  {
    return m_default;
  }

  /**
   * @brief check if param is set
   *
   * @return true
   * @return false
   */
  bool is_set()
  {
    return m_is_set;
  }

  /**
   * @brief get param value
   *
   * Support:
   * - Arithmetic types (int, float, double...)
   * - bool
   * - std::string
   *
   * Example:
   * @code
   * bc::param_t p;
   * int v = p->as<int>();
   * @endcode
   *
   * @endcode
   *
   * Extend for custom types:
   * @code
   * class MyClass
   * {
   * public:
   *   MyClass(std::string one, int i) {}
   * };
   *
   * template<>
   * MyClass bc::param::Param::as<MyClass>()
   * {
   *   return MyClass(m_str_value, 10);
   * }
   * @endcode
   * Then:
   * @code
   * bc::param_t p;
   * MyClass c = p->as<MyClass>();
   * @endcode
   *
   *
   * @tparam T
   * @return T
   */
  template<typename T>
  T as()
  {
    if constexpr(std::is_same_v<T, bool>)
      return m_is_set;
    else
      return utils::lexical_cast<T>(m_str_value);
  }

PRIVATE:
  Param(const std::string& name,
        const std::string& help,
        setter_fn_t setter,
        checker_fn_t checker)
    : m_raw_name(name), m_help(help), c_setter(setter)
  {
    for (auto& v : utils::split(name, '/'))
    {
      if (utils::is_short_param(v))
        m_short = utils::trim_param(v);
      else if (utils::is_long_param(v))
        m_long = utils::trim_param(v);
      else
      {
        ex::ExHandler::get().push(
          ex::InvalidParamError(
            name + " -> " + utils::wrap(v, "\"\"") + " is neither a valid short nor a valid long param")
         );
      }
      m_has_pname = true;
    }
    if (checker) c_checkers.push_back(checker);
  }

  const std::vector<std::tuple<checker_fn_t, param_t, checker_fn_t>>& get_dependency()
  {
    return m_depends_on;
  }

  const std::vector<std::tuple<checker_fn_t, param_t, check::checker_fn_t>>& get_banned()
  {
    return m_banned;
  }

  bool is_required()
  {
    return !m_has_default;
  }

  void set()
  {
    m_is_set = true;
    m_str_value = FLAG_VALUE;
  }

  std::string sp() const
  {
    if (m_short.empty())
      return m_short;
    return std::string("-") + m_short;
  }

  std::string lp() const
  {
    if (m_long.empty())
      return m_long;
    return std::string("--") + m_long;
  }

  bool hidden() const
  {
    return m_hidden;
  }

  Action get_action() const
  {
    return m_action;
  }

  const std::string& get_help() const
  {
    return m_help;
  }

  const std::string& get_meta() const
  {
    return m_meta;
  }

  bool has_short()
  {
    return !m_short.empty();
  }

  bool has_long()
  {
    return !m_long.empty();
  }

  bool has_both()
  {
    return has_short() && has_long();
  }

  bool is_flag()
  {
    return m_is_flag;
  }

  void process(const std::string& value)
  {
    m_str_value = value;
    if (!c_checkers.empty())
    {
      for (auto& cc : c_checkers)
      {
        auto [res, msg] = cc(m_raw_name, m_str_value);
        if (!res)
          throw ex::CheckFailedError(msg);
      }
      m_has_valid_value = true;
    }
    m_is_set = true;
    if (c_setter)
    {
      c_setter(value);
    }
    if (c_callback && m_callback_trigger && !m_as_default)
    {
      c_callback();
    }
  }

  void process_def()
  {
    m_as_default = true;
    process(m_str_value);
  }

PRIVATE:
  std::string m_raw_name  {};
  std::string m_str_value {};
  std::string m_help      {};
  std::string m_default   {};
  std::string m_short     {};
  std::string m_long      {};
  std::string m_meta      {conf::get().m_default_meta};

  Action m_action {Action::Nothing};

  std::vector<std::tuple<checker_fn_t, param_t, checker_fn_t>> m_depends_on {};
  std::vector<std::tuple<checker_fn_t, param_t, checker_fn_t>> m_banned {};

  bool m_has_default     {false};
  bool m_has_valid_value {false};
  bool m_has_pname       {false};
  bool m_is_set          {false};
  bool m_is_flag         {false};

  bool m_as_default      {false};
  bool m_hidden          {false};
  //CheckerMode m_check_mode {CheckerMode::AND};

PRIVATE:
  setter_fn_t      c_setter;
  checker_fn_t     c_checker;
  std::vector<checker_fn_t> c_checkers;
  callback_fn_t    c_callback;

  bool m_callback_trigger {false};
};

/**
 * @ingroup Param
 * @brief make a std::shared_ptr<Param>
 *
 * @param name parameter name -> "-p" || "--param" || "-p/--param"
 * @param help parameter help
 * @param setter a setter function
 * @param checker a checker function
 * @return param_t std::shared_ptr<Param>
 */
inline param_t make(const std::string& name,
                    const std::string& help,
                    setter_fn_t setter = nullptr,
                    checker_fn_t checker = nullptr)
{
  return std::make_shared<Param>(Param(name, help, setter, checker));
}

template<>
inline param_t Param::setter<setter_fn_t>(setter_fn_t& setter_callback)
{
  c_setter = setter_callback;
  return shared_from_this();
}

/**
 * @defgroup ParamGroup
 * @brief About bcli parameter groups
 */

class ParamGroup;

/**
 * @ingroup ParamGroup
 * @typedef pgroup_t
 * @brief pgroup_t
 *
 * ParamGroup constructor is private, ParamGroup is always used as std::shared_ptr<ParamGroup>
 */
using pgroup_t = std::shared_ptr<ParamGroup>;

/**
 * @ingroup ParamGroup
 * @brief ParamGroup
 *
 * ParamGroup class
 */
class ParamGroup
{
  friend pgroup_t make_group(const std::string& name,
                             const std::string& desc);
  friend class Command;
  friend class Parser<0>;
  friend class Parser<1>;

PRIVATE:
  ParamGroup(const std::string& name,
             const std::string& desc)
    : m_name(name), m_desc(desc)
  {}

  std::string get_help()
  {
    if (m_hidden)
      return "";

    std::stringstream help;
    const std::string& fsymb = conf::get().m_flag_symbol;
    help.fill(' ');
    help << utils::sp(2) << utils::wrap(m_name, "[]");
    if (!m_desc.empty()) help << " - " << m_desc;
    help << "\n";
    auto [maxs, maxl] = padding();
    for (auto& p : m_order)
    {
      if (p->hidden())
        continue;
      std::string fl = p->is_flag() ? " " + fsymb : "";
      help << utils::sp(4);
      if (!p->sp().empty())
        help << std::setw(maxs) << std::right << p->sp();
      else help << utils::sp(maxs);
      help << " ";
      if (!p->lp().empty())
        help << std::setw(maxl) << std::left << p->lp();
      else help << utils::sp(maxl-fl.size());

      help <<" - " << p->get_help() << " ";

      if (p->is_flag()) help << utils::wrap(fsymb, "[]");
      if (!p->get_def().empty()) help << utils::wrap(p->get_def(), "{}");
      help << "\n";
    }
    help << "\n";
    return help.str();
  }

public:
  /**
   * @brief add param
   *
   * @param name param name (ex: "-p/--param")
   * @param help parma help
   * @return param_t
   */
  param_t add_param(const std::string& name,
                    const std::string& help)
  {
    param_t p = make(name, help);
    add(p);
    return p;
  }

  /**
   * @brief get param_t
   *
   * @param pname param name
   * @return param_t
   */
  param_t get(const std::string& pname)
  {
    if (m_params.count(pname) > 0)
      return m_params.at(pname);
    else
      throw ex::UnknownParamError(pname + "doesn't exist in group " + m_name);
  }

  auto begin() { return m_order.begin(); }
  auto end() { return m_order.end(); }
  auto begin() const { return m_order.begin(); }
  auto end() const { return m_order.end(); }

  void hide()
  {
    m_hidden = true;
  }

  bool hidden()
  {
    return m_hidden;
  }

PRIVATE:
  bool exists(const std::string& pname)
  {
    std::string p = utils::trim_param(pname);
    return m_params.count(p) > 0;
  }

  const std::string& name()
  {
    return m_name;
  }

  bool empty()
  {
    return m_order.size() == 0;
  }

  void add(param_t p)
  {
    auto [sp, lp] = p->idx();
    already_exists(p, sp, lp);
    if (!sp.empty()) m_params.insert({sp, p});
    if (!lp.empty()) m_params.insert({lp, p});
    m_nbp++;
    m_order.push_back(p);
  }

  std::tuple<int, int> padding()
  {
    size_t max_s = 0, max_l = 0;
    for (auto& p: m_order)
    {
      max_s = p->sp().size() < max_s ? max_s : p->sp().size();
      max_l = p->lp().size() < max_l ? max_l : p->lp().size();
    }
    return std::make_tuple(max_s, max_l);
  }

  void already_exists(param_t p, const std::string& sp, const std::string& lp)
  {
    if (m_params.count(sp) > 0 || m_params.count(lp) > 0)
    {
      ex::ExHandler::get().push(
        ex::AlreadyExistsError(p->raw() + " already exists in group " + m_name));
    }
  }

PRIVATE:
  std::string                              m_name {};
  std::string                              m_desc {};
  std::unordered_map<std::string, param_t> m_params {};
  std::vector<param_t>                     m_order {};
  size_t                                   m_nbp {0};
  bool                                     m_hidden {false};
};

/**
 * @ingroup ParamGroup
 * @brief make a std::shared_ptr<ParamGroup>
 *
 * @param name group name
 * @param desc group description
 * @return pgroup_t
 */
inline pgroup_t make_group(const std::string& name,
                           const std::string& desc)
{
  return std::make_shared<ParamGroup>(ParamGroup(name, desc));
}

/**
 * @defgroup Command
 * @brief About bcli commands
 */

class Command;

/**
 * @ingroup Command
 * @typedef cmd_t
 * @brief Command constuctor is private, Command is always used as std::shared_ptr<Command>
 *
 */
using cmd_t = std::shared_ptr<Command>;

class Commands;
/**
 * @ingroup Command
 * @brief Command
 *
 * Command class
 */
class Command
{
  friend cmd_t make_cmd(const std::string& name,
                        const std::string& desc,
                        help_fn_t help);

  friend class Parser<0>;
  friend class Parser<1>;
  friend class Commands;

PRIVATE:
  Command(const std::string& name,
          const std::string& desc,
          help_fn_t help = nullptr)
    : m_name(name), m_desc(desc), m_help(help)
  {}

public:
  /**
   * @brief set command help
   *
   * @param help
   */
  void set_help(help_fn_t help)
  {
    m_help = help;
  }

  /**
   * @brief add param to current group
   *
   * @param name parameter name (ex: "-p/--param")
   * @param help  parameter help
   * @return param_t
   */
  param_t add_param(const std::string& name,
                   const std::string& help)
  {
    return m_current_group->add_param(name, help);
  }

  /**
   * @brief add group to command
   *
   * @param name group name
   * @param desc group description
   * @return pgroup_t
   */
  pgroup_t add_group(const std::string& name,
                     const std::string& desc)
  {
    param::pgroup_t group = make_group(name, desc);
    add(group);
    return group;
  }

  /**
   * @brief add a group with common parameters
   *
   * "-h/--help" -> trigger help \n
   * "-v/--verbose" -> trigger version \n
   * "-d/--debug" \n
   * "--version" \n
   * @see bc::Action
   * @see Configuration
   *
   * @param name group name
   * @return pgroup_t
   */
  pgroup_t add_common(const std::string& name = "common")
  {
    if (conf::get().has_common())
    {
      pgroup_t grp = add_group(name, {});
      if (conf::get().m_help)
        grp->add_param("-h/--help", "Show this message and exit.")
           ->as_flag()->action(Action::ShowHelp);
      if (conf::get().m_verbose)
        grp->add_param("-v/--verbose", "Verbose mode.")->as_flag();
      if (conf::get().m_debug)
        grp->add_param("-d/--debug", "Debug mode.")->as_flag();
      if (conf::get().m_version)
       grp->add_param("--version", "Show version and exit.")
           ->as_flag()->action(Action::ShowVersion);
      return grp;
    }
    else
      return nullptr;
  }

  /**
   * @brief get parameter
   *
   * @code
   * cmd->add_param("-p/--param", "help param");
   *
   * cmd->getp("p");
   * cmd->getp("param");
   * cmd->getp("-p");
   * cmd->getp("--param");
   *
   * @endcode
   * @param name
   * @return param_t
   */
  param_t getp(const std::string& name)
  {
    std::string name_ = utils::trim_param(name);
    for (auto& grp : m_order)
      if (grp->exists(name_))
        return grp->get(name_);
    return nullptr;
  }

  /**
   * @brief get positional parameters
   *
   * @return const std::vector<std::string>&
   */
  const std::vector<std::string>& get_positionals() const
  {
    return m_positionals;
  }

  /**
   * @brief set positionals help
   *
   * @param usage
   * @param help
   */
  void set_positionals_help(const std::string& usage, const std::string& help)
  {
    m_usage_pos = usage;
    m_help_pos = help;
  }

  /**
   * @brief set number of positionals
   *
   * @param nb number of positionals
   * @param usage
   * @param help
   */
  void set_positionals(size_t nb,
                       const std::string& usage,
                       const std::string& help)
  {
    m_e_pos = nb;
    m_usage_pos = usage;
    m_help_pos = help;
    m_checkp = true;
  }

  /**
   * @brief set positionals bounds
   *
   * @param min min number of positionals
   * @param max max number of positionals
   * @param usage
   * @param help
   */
  void set_positional_bounds(size_t min, size_t max,
                             const std::string& usage,
                             const std::string& help)
  {
    m_l_pos = min;
    m_u_pos = max;
    m_usage_pos = usage;
    m_help_pos = help;
    m_checkp = true;
    m_bmode = true;
  }

  /**
   * @brief set checker for positionals
   *
   * @param checker
   */
  void positionals_checker(checker_fn_t checker)
  {
    c_pchecker = checker;
  }

  void positionals_setter(setter_fn_t setter)
  {
    c_psetter = setter;
  }

  auto begin() { return m_order.begin(); }
  auto end() { return m_order.end(); }
  auto begin() const { return m_order.begin(); }
  auto end() const { return m_order.end(); }

PRIVATE:
  std::string& name()
  {
    return m_name;
  }

  std::string& desc()
  {
    return m_desc;
  }

  std::string get_help(const std::string& main_name, const std::string& version, bool cmd_mode)
  {
    if (m_help)
      return m_help();
    else
    {
      std::stringstream ss;
      ss << main_name << " ";
      if (cmd_mode) ss << m_name << " ";
      ss << version << "\n\n";

      ss << "DESCRIPTION" << "\n" << utils::sp(2);
      ss << m_desc << "\n\n";

      ss << "USAGE" << "\n" << utils::sp(2);
      if (cmd_mode) ss << main_name << " ";
      ss << m_name << " " << get_usage(main_name, cmd_mode) << "\n\n";
      if (!m_help_pos.empty())
      {
        ss << "POSITIONALS" << "\n";
        ss << utils::sp(2) << m_help_pos << "\n\n";
      }
      ss << "OPTIONS" << "\n";
      for (auto& g : m_order)
        if (!g->empty())
          ss << g->get_help();
      return ss.str();
    }
  }

  std::string get_usage(const std::string& main_name, bool cmd_mode)
  {
    std::stringstream ss;
    std::vector<std::string> required;
    std::vector<std::string> optionals;
    std::vector<std::string> flags;

    for (auto& group : m_order)
    {
      if (group->hidden())
        continue;
      for (auto& p : *group)
      {
        if (p->hidden())
          continue;
        std::string bds = p->is_required() ? "" : "[]";
        std::string raw;
        if (p->is_flag())
          raw = p->raw();
        else
          raw = p->raw() + " " + utils::wrap(p->get_meta(), "<>");
        std::string fm = utils::wrap(raw, bds);

        if (p->is_required())
          required.push_back(fm);
        else if (p->is_flag())
          flags.push_back(fm);
        else
          optionals.push_back(fm);
      }
    }

    size_t pad;
    if (cmd_mode) pad = main_name.size() + 1 + m_name.size();
    else pad = m_name.size();

    auto format = [&ss, &pad](const std::vector<std::string> vec, int c) -> int {
      size_t max_size = 90;
      int current = c;
      for (auto& f : vec)
      {
        if (current + f.size() < max_size)
        {
          ss << f << " ";
          current += f.size();
        }
        else
        {
          ss << "\n" << utils::sp(pad+3) << f << " ";
          current = f.size()+pad+3;
        }
      }
      return current;
    };
    format(flags, format(optionals, format(required, 0)));
    ss << m_usage_pos;
    return ss.str();
  }

  std::tuple<bool, std::string> check_positionals()
  {
    if (m_checkp)
    {
      if (m_bmode)
      {
        if (m_positionals.size() < m_l_pos)
          return std::make_tuple(
            false,
            "requires at least " + std::to_string(m_l_pos) + " positionals.");
        else if (m_positionals.size() > m_u_pos)
          return std::make_tuple(
            false,
            "requires at most " + std::to_string(m_u_pos) + " positionals.");
      }
      else
        return std::make_tuple(
          m_positionals.size() == m_e_pos,
          "number of positionals must be " + std::to_string(m_e_pos)
        );
    }

    if (c_pchecker)
    {
      int i = 0;
      for (auto& p : m_positionals)
      {
        auto [res, msg] = c_pchecker("positionals" + utils::wrap(std::to_string(i), "[]"), p);
        if (!res)
          return std::make_tuple(res, msg);
        i++;
      }
    }

    return std::make_tuple(true, "");
  }

  void push_positionals(const std::string& arg)
  {
    m_positionals.push_back(arg);
    m_nb_pos++;
    if (c_psetter)
      c_psetter(arg);
  }

  void add(pgroup_t pg)
  {
    if (m_groups.count(pg->name()) > 0)
    {
      ex::ExHandler::get().push(
        ex::AlreadyExistsError(pg->name() + "group already exists in command " + m_name));
    }
    else
    {
      m_groups.insert({pg->name(), pg});
      m_order.push_back(pg);
      m_current_group = pg;
    }
  }

  pgroup_t get(const std::string& name)
  {
    if (m_groups.count(name) == 0)
      throw ex::UnknownGroupError(name + "doesn't exist in command " + m_name);
    else
      return m_groups.at(name);
  }

PRIVATE:
  std::string m_name {};
  std::string m_desc {};
  help_fn_t   m_help {nullptr};
  pgroup_t    m_def_group {make_group("", "")};
  pgroup_t    m_current_group {};
  bool        m_def_use {false};
  size_t      m_nb_pos {0};
  std::unordered_map<std::string, pgroup_t> m_groups;
  std::vector<pgroup_t>                     m_order;
  std::vector<std::string> m_positionals;

  std::string m_help_pos {};
  std::string m_usage_pos {};
  size_t m_e_pos {0};
  size_t m_l_pos {0};
  size_t m_u_pos {0};
  bool   m_checkp {false};
  bool   m_bmode {false};

  checker_fn_t c_pchecker {nullptr};
  setter_fn_t  c_psetter {nullptr};
};

/**
 * @ingroup Command
 * @brief make a std::shared_ptr<Command>
 *
 * @param name command name
 * @param desc command description
 * @param help command help
 * @return cmd_t
 */
inline cmd_t make_cmd(const std::string& name,
                      const std::string& desc,
                      help_fn_t help = nullptr)
{
  return std::make_shared<Command>(Command(name, desc, help));
}

class Commands;
/**
 * @ingroup Command
 * @typedef cmds_t
 * @brief Commands constuctor is private, Commands is always used as std::shared_ptr<Commands>
 *
 */
using cmds_t = std::unique_ptr<Commands>;

class Commands
{
  friend cmds_t make_cmds(const std::string&, const std::string&, const std::string&);
  friend class Parser<1>;
  friend class Parser<0>;

PRIVATE:
  Commands(const std::string& name, const std::string& desc, const std::string& version)
    : m_name(name), m_desc(desc), m_version(version)
  {}

  void set_help(help_fn_t help)
  {
    m_help = help;
  }

  cmd_t get(const std::string& name)
  {
    return m_cmds.at(name);
  }

  bool exists(const std::string& name)
  {
    return m_cmds.count(name) > 0;
  }

  std::vector<std::string> list()
  {
    std::vector<std::string> ret;
    for (auto& c : m_order)
      ret.push_back(c->name());
    return ret;
  }

  std::string get_help()
  {
    if (m_help)
      return m_help();
    else
    {
      size_t pad = (*std::max_element(m_order.begin(), m_order.end(), [](cmd_t a, cmd_t b) -> bool {
          return a->name().size() < b->name().size();
      }))->name().size();
      std::stringstream help;
      help << m_name << " " << m_version << "\n\n";
      help << "DESCRIPTION" << "\n" << utils::sp(2) << m_desc << "\n\n";
      help << "USAGE" << "\n" << utils::sp(2);
      help << m_name << " " << utils::wrap(utils::join(list(), "|"), "[]") << "\n\n";
      help << "COMMANDS" << "\n";
      for (auto& cmd: m_order)
      {
        help << utils::sp(2) << std::setw(pad) << std::left << cmd->name() << " - " << cmd->desc() << "\n";
      }
      return help.str();
    }
  }

  void add(cmd_t c)
  {
    m_cmds.insert({c->name(), c});
    m_order.push_back(c);
  }

PRIVATE:
  std::string m_name;
  std::string m_desc;
  std::string m_version;
  std::unordered_map<std::string, cmd_t> m_cmds;
  std::vector<cmd_t>                     m_order;
  help_fn_t                              m_help {nullptr};
};

inline cmds_t make_cmds(const std::string& name,
                        const std::string& desc,
                        const std::string& version)
{
  return std::make_unique<Commands>(Commands(name, desc, version));
}

} // end of namespace param

/**
 * @defgroup Parser
 * @brief About bcli parser.
 *
 */

using pgroup_t = param::pgroup_t;
using param_t  = param::param_t;
using cmd_t    = param::cmd_t;
using conf     = config::Config;

#define ENABLE_IF(n)                                               \
template<int M = Mode,                                             \
         typename = typename std::enable_if<M == n, void>::type>   \

#define ENABLE_IF_P(n, pack)                                       \
template<pack,                                                     \
         int M = Mode,                                             \
         typename = typename std::enable_if<M == n, void>::type>   \

#define BCLI_PARSE(parser, argc, argv) \
try {                                  \
  parser.parse(argc, argv);            \
} catch (const bc::ex::BCliError& e) { \
  bc::utils::exit_bcli(e);             \
  exit(EXIT_FAILURE);                  \
}                                      \

/**
 * @ingroup Parser
 * @brief bcli parser
 *
 * @tparam Mode 0: classical, 1: multiple commands
 */
template<int Mode = 0>
class Parser
{
public:
  Parser(const std::string& name,
         const std::string& desc,
         const std::string& version,
         const std::string& authors = {})
    : m_name(name), m_desc(desc), m_version(version), m_authors(authors)
  {
    m_current_cmd->add(m_current_group);
  }

  /**
   * @ingroup Parser
   * @brief parse argv
   *
   * @param argc
   * @param argv
   */
  void parse(int argc, char* argv[])
  {
    ex::ExHandler::get().check();

    if (!m_is_cmd_mode || m_bypass)
    {
      for (int i=1; i<argc; i++)
      {
        switch (process_arg(argv[i]))
        {
        case Action::ShowHelp:
          show_help();
          throw ex::BCliError("", "", ex::ExitCodes::Failure);
        case Action::ShowVersion:
          show_version();
          throw ex::BCliError("", "", ex::ExitCodes::Failure);
        default: break;
        }
      }
      if (m_is_param)
        throw ex::MissingValueError(m_current + " needs a value.");
      check_consistency();
    }
    else
    {
      if (argc < 2)
      {
        show_help();
        throw ex::BCliError("", "", ex::ExitCodes::Failure);
      }
      if (m_cmds->exists(argv[1]))
      {
        m_current_cmd = m_cmds->get(argv[1]);
        m_bypass = true;
        parse(argc-1, argv+1);
      }
      else
      {
        std::stringstream ss;
        ss << "Unknown command: " << argv[1] << ", ";
        ss << "choices -> "<< utils::wrap(utils::join(m_cmds->list(), "|"), "[]");
        throw ex::UnknownCmdError(ss.str());
      }
    }
  }

public:
  /**
   * @ingroup Parser
   * @brief add commmand
   *
   * @param name command name
   * @param desc command description
   * @param help command help
   * @return param::cmd_t
   */
  ENABLE_IF(1)
  param::cmd_t add_command(const std::string& name,
                     const std::string& desc,
                     param::help_fn_t help = nullptr)
  {
    param::cmd_t cmd = param::make_cmd(name, desc, help);
    m_current_cmd = cmd;
    m_current_cmd->add(param::make_group(conf::get().m_default_grp, ""));
    m_cmds->add(cmd);
    m_is_cmd_mode = true;
    return cmd;
  }

  ENABLE_IF(1)
  bool is(const std::string& cmd_name)
  {
    return (cmd_name == m_current_cmd->name());
  }

  /**
   * @ingroup Parser
   * @brief add group
   *
   * @param name group name
   * @param desc group description
   * @return param::pgroup_t
   */
  ENABLE_IF(0)
  param::pgroup_t add_group(const std::string& name,
                            const std::string& desc)
  {
    if (m_current_group != m_current_cmd->get(m_current_group->name()))
      m_current_cmd->add(m_current_group);
    m_current_group = param::make_group(name, desc);
    m_current_cmd->add(m_current_group);
    return m_current_group;
  }

  /**
   * @ingroup Parser
   * @brief add common
   *
   * @param name name
   * @return param::pgroup_t
   */
  ENABLE_IF(0)
  param::pgroup_t add_common(const std::string& name = "common")
  {
    m_current_cmd->add_common(name);
    m_current_group = m_current_cmd->get(name);
    return m_current_group;
  }

  /**
   * @ingroup Parser
   * @brief add param
   *
   * @param name param name
   * @param help param help
   * @return param::param_t
   */
  ENABLE_IF(0)
  param::param_t add_param(const std::string& name,
                           const std::string& help)
  {
    param::param_t p = param::make(name, help);
    m_current_group->add(p);
    return p;
  }

  /**
   * @ingroup Parser
   * @brief get positional arguments
   *
   * @return const std::vector<std::string>&
   */
  const std::vector<std::string>& get_positionals() const
  {
    return m_current_cmd->get_positionals();
  }

  /**
   * @ingroup Parser
   * @brief set number of positionals
   *
   * @param nb number of positionals
   * @param usage positionals usage
   * @param help positionals help
   */
  ENABLE_IF(0)
  void set_positionals(size_t nb, const std::string& usage, const std::string& help)
  {
    m_current_cmd->set_positionals(nb, usage, help);
  }

  /**
   * @ingroup Parser
   * @brief set positional bounds
   *
   * @param min min number of positionals
   * @param max max number of positionals
   * @param usage positionals usage
   * @param help positionals help
   */
  ENABLE_IF(0)
  void set_positional_bounds(size_t min, size_t max,
                             const std::string& usage,
                             const std::string& help)
  {
    m_current_cmd->set_positional_bounds(min, max, usage, help);
  }

  /**
   * @ingroup Parser
   * @brief set positionals usage and help
   *
   * @param usage positionals usage
   * @param help positionals help
   */
  ENABLE_IF(0)
  void set_positionals_help(const std::string& usage, const std::string& help)
  {
    m_current_cmd->set_positionals_help(usage, help);
  }

  /**
   * @ingroup Parser
   * @brief set positionals checker
   *
   * @param checker
   */
  ENABLE_IF(0)
  void positionals_checker(check::checker_fn_t checker)
  {
    m_current_cmd->positionals_checker(checker);
  }

  /**
   * @ingroup Parser
   * @brief set help
   *
   * @param help
   */
  void set_help(std::function<std::string()> help)
  {
    if constexpr(Mode == 0)
      m_current_cmd->set_help(help);
    else
      m_cmds->set_help(help);
  }

  /**
   * @ingroup Parser
   * @brief get param
   *
   * For "-p/--param":
   * @code
   * cli->getp("p");
   * cli->getp("-p");
   * cli->getp("param");
   * cli->getp("--param");
   * @endcode
   * @param pname param name
   * @return param_t
   */
  param_t getp(const std::string& pname) const
  {
    std::string pname_ = utils::trim_param(pname);
    for (auto& group : *m_current_cmd)
    {
      if (group->exists(pname_))
        return group->get(pname_);
    }
    return nullptr;
  }

  void show_help()
  {
    if (!m_is_cmd_mode || m_bypass)
    {
      std::cerr << m_current_cmd->get_help(m_name, m_version, m_bypass) << std::endl;
    }
    else
    {
      std::cerr << m_cmds->get_help() << std::endl;
    }
  }

PRIVATE:
  void show_version()
  {
    std::cerr << m_name << " " << m_version << std::endl;
  }

  bool pexists(const std::string& pname) const
  {
    for (auto& group : *m_current_cmd)
    {
      if (group->exists(pname))
        return true;
    }
    return false;
  }

  Action process_arg(const std::string& arg)
  {
    if (utils::is_param(arg))
    {
      if (!pexists(arg))
        throw ex::InvalidParamError("Unknown param: " + arg + ".");
      else if (m_is_param)
        throw ex::MissingValueError(m_current + "needs a value.");
      m_current = arg;
      m_is_param = true;

      param::param_t cp = get_current_param(m_current);
      if (cp->is_flag())
      {
        m_is_param = false;
        m_last_is_flag = true;
        cp->set();
        get_current_param(m_current)->process(FLAG_VALUE);
      }

      if (cp->get_action() != Action::Nothing) return cp->get_action();
    }
    else if (!m_is_param)
    {
      m_current_cmd->push_positionals(arg);
    }
    else if (m_is_param)
    {
      std::string parg;
      if (utils::startswith(arg, "[-") && utils::endswith(arg, "]"))
        parg = utils::unwrap(arg, 1);
      else
        parg = arg;
      get_current_param(m_current)->process(parg);
      m_is_param = false;
    }
    return Action::Nothing;
  }

  param::param_t get_current_param(const std::string& p) const
  {
    std::string pname = utils::trim_param(p);
    for (auto& group: *m_current_cmd)
      if (group->exists(pname))
        return group->get(pname);
    return nullptr;
  }

  void check_consistency()
  {
    for (auto& group: *m_current_cmd)
    {
      for (auto& p : *group)
      {
        if (p->is_required() && !p->is_set())
          throw ex::RequiredParamError(p->raw() + " is required.");
        else if (!p->is_flag() && !p->get_def().empty())
          p->process_def();

        for (auto& [c, d, dc] : p->get_dependency())
        {
          auto [res, msg] = c(p->raw(), p->value());
          if (!dc)
          {
            if (res)
              if (!d->is_set())
                throw ex::DependsError(utils::format_depend_errors(p, d, msg));
          }
          else
          {
            auto [dres, dmsg] = dc(d->raw(), d->value());
            if (res)
              if (!dres)
                throw ex::DependsError(utils::format_depend_errors(p, d, dmsg));
          }
        }

        for (auto& [c, d, dc] : p->get_banned())
        {
          auto [res, msg] = c(p->raw(), p->value());
          if (!dc)
          {
            if (res)
              if (d->is_set())
                throw ex::BannedError(utils::format_banned_errors(p, d, msg));
          }
          else
          {
            auto [dres, dmsg] = dc(d->raw(), d->value());
            if (res)
              if (dres)
                throw ex::BannedError(utils::format_banned_errors(p, d, dmsg));
          }
        }
      }
    }
    auto [res, msg] = m_current_cmd->check_positionals();
    if (!res)
      throw ex::PositionalsError(msg);
    return;
  }

PRIVATE:
  std::string m_name {};
  std::string m_desc {};
  std::string m_version {};
  std::string m_authors {};

  param::pgroup_t m_current_group {param::make_group(conf::get().m_default_grp,
                                                     "global parameters")};
  param::cmd_t    m_current_cmd {param::make_cmd(m_name, m_desc)};
  param::cmds_t   m_cmds {param::make_cmds(m_name, m_desc, m_version)};

  std::string m_current {};

  std::vector<ex::BCliError> m_impl_exceptions {};
  std::vector<ex::BCliError> m_usage_exceptions {};

  bool m_is_param {false};
  bool m_last_is_flag {false};
  bool m_is_cmd_mode {false};
  bool m_bypass {false};
};

} // end of namespace bc

#undef PRIVATE
