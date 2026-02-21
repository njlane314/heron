/* -- C++ -- */
/**
 *  @file io/src/Dataset.cpp
 *
 *  @brief Implementation for Dataset sample-list loading.
 */

#include "Dataset.hh"

#include <fstream>
#include <stdexcept>
#include <sstream>

namespace heron
{

void Dataset::load_samples(const std::string &sample_list_path)
{
    this->clear();

    std::ifstream fin(sample_list_path);
    if (!fin)
    {
        throw std::runtime_error(std::string("Dataset::load_samples: failed to open sample list: ") + sample_list_path);
    }

    std::string line;
    while (std::getline(fin, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        std::istringstream row(line);
        std::string sample_name;
        std::string sample_origin;
        std::string beam_mode;
        std::string output_path;

        if (!std::getline(row, sample_name, '\t') ||
            !std::getline(row, sample_origin, '\t') ||
            !std::getline(row, beam_mode, '\t') ||
            !std::getline(row, output_path))
        {
            continue;
        }

        if (sample_name == "sample_name" || output_path.empty())
        {
            continue;
        }

        m_samples.push_back(SampleIO::read(output_path));
    }

}

void Dataset::clear()
{
    m_samples.clear();
}

const std::vector<Dataset::Sample> &Dataset::samples() const
{
    return m_samples;
}

} // namespace heron
