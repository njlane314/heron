/* -- C++ -- */
/**
 *  @file io/include/Dataset.hh
 *
 *  @brief Dataset container that loads sample metadata from a sample list.
 */

#ifndef heron_io_dataset_h
#define heron_io_dataset_h

#include <string>
#include <vector>

#include "SampleIO.hh"

namespace heron
{

/** \brief Container for loaded analysis samples. */
class Dataset
{
  public:
    using Sample = SampleIO::Sample;

    /// Load samples from a TSV sample list.
    bool load_samples(const std::string &sample_list_path);

    /// Clear loaded samples.
    void clear();

    /// Access loaded samples.
    const std::vector<Sample> &samples() const;

  private:
    std::vector<Sample> m_samples;
};

} // namespace heron

#endif // heron_io_dataset_h
