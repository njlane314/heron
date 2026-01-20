/* -- C++ -- */
/**
 *  @file  io/include/SampleRootIO.hh
 *
 *  @brief ROOT IO helpers for aggregated samples.
 */

#ifndef NUXSEC_SAMPLE_SAMPLE_ROOT_IO_H
#define NUXSEC_SAMPLE_SAMPLE_ROOT_IO_H

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

#endif // NUXSEC_SAMPLE_SAMPLE_ROOT_IO_H
