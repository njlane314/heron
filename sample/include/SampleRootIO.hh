/* -- C++ -- */
/**
 *  @file  sample/include/SampleRootIO.hh
 *
 *  @brief ROOT IO helpers for aggregated sample records.
 */

#ifndef Nuxsec_SAMPLE_SAMPLE_ROOT_IO_H_INCLUDED
#define Nuxsec_SAMPLE_SAMPLE_ROOT_IO_H_INCLUDED

#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>

#include <string>

#include "Sample.hh"

namespace nuxsec
{

class SampleRootIO
{
  public:
    static void write(const Sample &sample, const std::string &out_file);
    static Sample read(const std::string &in_file);

  private:
    static constexpr const char *kSampleDir = "nuxsec_sample";
};

} // namespace nuxsec

#endif // Nuxsec_SAMPLE_SAMPLE_ROOT_IO_H_INCLUDED
