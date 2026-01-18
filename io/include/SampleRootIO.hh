/* -- C++ -- */
/**
 *  @file  io/include/SampleRootIO.hh
 *
 *  @brief ROOT IO helpers for aggregated samples.
 */

#ifndef Nuxsec_IO_SAMPLE_ROOT_IO_H_INCLUDED
#define Nuxsec_IO_SAMPLE_ROOT_IO_H_INCLUDED

#include <string>

#include "Sample.hh"

namespace nuxsec
{

class SampleRootIO
{
  public:
    static void write(const Sample &sample, const std::string &out_file);
    static Sample read(const std::string &in_file);
};

} // namespace nuxsec

#endif // Nuxsec_IO_SAMPLE_ROOT_IO_H_INCLUDED
