/* -- C++ -- */
/**
 *  @file  framework/core/src/Dataset.cc
 *
 *  @brief Backend dataset loader implementation for event workflows.
 */

#include "Dataset.hh"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <utility>

std::vector<SampleListEntry> Dataset::read_samples(const std::string &list_path,
                                                   const bool allow_missing,
                                                   const bool require_nonempty)
{
    std::ifstream fin(list_path);
    if (!fin)
    {
        if (allow_missing && errno == ENOENT)
        {
            return {};
        }
        throw std::runtime_error("Failed to open sample list: " + list_path +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }

    std::vector<SampleListEntry> entries;
    std::string line;
    bool first_nonempty = true;
    while (std::getline(fin, line))
    {
        line = trim(line);
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        const auto fields = split_tabs(line);
        if (fields.size() < 4)
        {
            throw std::runtime_error("Malformed sample list entry: " + line);
        }

        if (first_nonempty && fields[0] == "sample_name")
        {
            first_nonempty = false;
            continue;
        }
        first_nonempty = false;

        SampleListEntry entry;
        entry.sample_name = fields[0];
        entry.sample_origin = fields[1];
        entry.beam_mode = fields[2];
        entry.output_path = fields[3];
        entries.push_back(std::move(entry));
    }

    if (require_nonempty && entries.empty())
    {
        throw std::runtime_error("Sample list is empty: " + list_path);
    }

    return entries;
}

Dataset Dataset::load(const std::string &list_path)
{
    const auto entries = read_samples(list_path);

    Dataset dataset;
    dataset.m_inputs.reserve(entries.size());
    dataset.m_sample_infos.reserve(entries.size());

    for (const auto &entry : entries)
    {
        SampleIO::Sample sample = SampleIO::read(entry.output_path);

        nu::SampleInfo info;
        info.sample_name = sample.sample_name;
        info.sample_rootio_path = entry.output_path;
        info.sample_origin = static_cast<int>(sample.origin);
        info.beam_mode = static_cast<int>(sample.beam);
        info.subrun_pot_sum = sample.subrun_pot_sum;
        info.db_tortgt_pot_sum = sample.db_tortgt_pot_sum;
        info.db_tor101_pot_sum = sample.db_tor101_pot_sum;
        dataset.m_sample_infos.push_back(std::move(info));

        DatasetInput input;
        input.entry = entry;
        input.sample = std::move(sample);
        dataset.m_inputs.push_back(std::move(input));
    }

    return dataset;
}
