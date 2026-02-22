/* -- C++ -- */
/**
 *  @file  framework/core/include/Dataset.hh
 *
 *  @brief Backend dataset loader for event workflows, handling sample IO and
 *         metadata mapping from sample lists.
 */
#ifndef HERON_CORE_DATASET_H
#define HERON_CORE_DATASET_H

#include <string>
#include <vector>

#include "EventListIO.hh"
#include "SampleCLI.hh"
#include "SampleIO.hh"

struct DatasetInput
{
    SampleListEntry entry;
    SampleIO::Sample sample;
};

class Dataset
{
  public:
    static Dataset load(const std::string &list_path);

    const std::vector<DatasetInput> &inputs() const noexcept { return m_inputs; }
    const std::vector<nu::SampleInfo> &sample_infos() const noexcept { return m_sample_infos; }

  private:
    static std::vector<SampleListEntry> read_samples(const std::string &list_path,
                                                     bool allow_missing = false,
                                                     bool require_nonempty = true);

    std::vector<DatasetInput> m_inputs;
    std::vector<nu::SampleInfo> m_sample_infos;
};

#endif // HERON_CORE_DATASET_H
