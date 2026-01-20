/* -- C++ -- */
/**
 *  @file  io/include/TemplateRootIO.hh
 *
 *  @brief ROOT IO helpers for template histogram output.
 */

#ifndef NUXSEC_IO_TEMPLATE_ROOT_IO_H
#define NUXSEC_IO_TEMPLATE_ROOT_IO_H

#include <string>
#include <utility>
#include <vector>

class TH1;

namespace nuxsec
{

struct TemplateWriteOptions
{
    std::string top_dir = "nuxsec/results";
    bool overwrite = true;
};

class TemplateRootIO
{
  public:
    static void write_histograms(const std::string &root_path,
                                 const std::string &sample_name,
                                 const std::vector<std::pair<std::string, const TH1 *>> &hists,
                                 const TemplateWriteOptions &opt = {});

    static void write_string_meta(const std::string &root_path,
                                  const std::string &sample_name,
                                  const std::string &key,
                                  const std::string &value,
                                  const TemplateWriteOptions &opt = {});

    static void write_double_meta(const std::string &root_path,
                                  const std::string &sample_name,
                                  const std::string &key,
                                  double value,
                                  const TemplateWriteOptions &opt = {});

    static void write_syst_histograms(const std::string &root_path,
                                      const std::string &sample_name,
                                      const std::string &syst_name,
                                      const std::string &variation,
                                      const std::vector<std::pair<std::string, const TH1 *>> &hists,
                                      const TemplateWriteOptions &opt = {});

    static void write_syst_flag_meta(const std::string &root_path,
                                     const std::string &syst_name,
                                     const std::string &key,
                                     const std::string &value,
                                     const TemplateWriteOptions &opt = {});
};

} // namespace nuxsec

#endif // NUXSEC_IO_TEMPLATE_ROOT_IO_H
