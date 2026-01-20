/* -- C++ -- */
#ifndef Nuxsec_APPS_APPSAMPLECOMMAND_H_INCLUDED
#define Nuxsec_APPS_APPSAMPLECOMMAND_H_INCLUDED

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "AppUtils.hh"
#include "SampleAggregator.hh"
#include "SampleRootIO.hh"
#include "SampleTypes.hh"

namespace nuxsec
{

namespace app
{

struct SampleArgs
{
    std::string sample_name;
    std::string filelist_path;
    std::string output_path;
    std::string sample_list_path;
};

inline SampleArgs parse_sample_spec(const std::string &spec)
{
    const auto pos = spec.find(':');
    if (pos == std::string::npos)
    {
        throw std::runtime_error("Bad sample spec (expected NAME:FILELIST): " + spec);
    }

    SampleArgs out;
    out.sample_name = nuxsec::app::trim(spec.substr(0, pos));
    out.filelist_path = nuxsec::app::trim(spec.substr(pos + 1));

    if (out.sample_name.empty() || out.filelist_path.empty())
    {
        throw std::runtime_error("Bad sample spec: " + spec);
    }

    out.output_path = "build/out/sample/sample_root_" + out.sample_name + ".root";
    out.sample_list_path = "build/out/sample/samples.tsv";

    return out;
}

inline SampleArgs parse_sample_args(const std::vector<std::string> &args, const std::string &usage)
{
    if (args.size() != 1)
    {
        throw std::runtime_error(usage);
    }

    return parse_sample_spec(args[0]);
}

inline void update_sample_list(const std::string &list_path,
                               const nuxsec::Sample &sample,
                               const std::string &output_path)
{
    auto entries = nuxsec::app::read_sample_list(list_path, true, false);
    const std::string kind_name = nuxsec::sample_kind_name(sample.kind);
    const std::string beam_name = nuxsec::beam_mode_name(sample.beam);

    bool updated = false;
    for (auto &entry : entries)
    {
        if (entry.sample_name == sample.sample_name &&
            entry.sample_kind == kind_name &&
            entry.beam_mode == beam_name)
        {
            entry.output_path = output_path;
            updated = true;
            break;
        }
    }
    if (!updated)
    {
        nuxsec::app::SampleListEntry entry;
        entry.sample_name = sample.sample_name;
        entry.sample_kind = kind_name;
        entry.beam_mode = beam_name;
        entry.output_path = output_path;
        entries.push_back(std::move(entry));
    }

    nuxsec::app::write_sample_list(list_path, std::move(entries));
}

inline int run_sample(const SampleArgs &sample_args, const std::string &log_prefix)
{
    const std::string db_path = "/exp/uboone/data/uboonebeam/beamdb/run.db";
    const auto files = nuxsec::app::read_file_list(sample_args.filelist_path);

    std::filesystem::path output_path(sample_args.output_path);
    if (!output_path.parent_path().empty())
    {
        std::filesystem::create_directories(output_path.parent_path());
    }
    std::filesystem::path sample_list_path(sample_args.sample_list_path);
    if (!sample_list_path.parent_path().empty())
    {
        std::filesystem::create_directories(sample_list_path.parent_path());
    }

    nuxsec::Sample sample = nuxsec::SampleAggregator::aggregate(sample_args.sample_name, files, db_path);
    nuxsec::SampleRootIO::write(sample, sample_args.output_path);
    update_sample_list(sample_args.sample_list_path, sample, sample_args.output_path);

    std::cerr << "[" << log_prefix << "] sample=" << sample.sample_name
              << " fragments=" << sample.fragments.size()
              << " pot_sum=" << sample.subrun_pot_sum
              << " db_tortgt_pot_sum=" << sample.db_tortgt_pot_sum
              << " normalisation=" << sample.normalisation
              << " normalised_pot_sum=" << sample.normalised_pot_sum
              << " output=" << sample_args.output_path
              << " sample_list=" << sample_args.sample_list_path
              << "\n";

    return 0;
}

}

}

#endif
