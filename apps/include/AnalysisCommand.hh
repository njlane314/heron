/* -- C++ -- */
#ifndef NUXSEC_APPS_ANALYSISCOMMAND_H
#define NUXSEC_APPS_ANALYSISCOMMAND_H

#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TH1.h>
#include <TH1D.h>

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>

#include "AnalysisConfigService.hh"
#include "AnalysisIO.hh"
#include "ColumnDerivationService.hh"
#include "RDataFrameFactory.hh"
#include "SampleIO.hh"
#include "Utils.hh"

namespace nuxsec
{

namespace app
{

struct AnalysisArgs
{
    std::string list_path;
    std::string output_root;
    int nthreads = 1;
};

inline AnalysisArgs parse_analysis_args(const std::vector<std::string> &args, const std::string &usage)
{
    if (args.size() < 2 || args.size() > 3)
    {
        throw std::runtime_error(usage);
    }

    AnalysisArgs out;
    out.list_path = nuxsec::app::trim(args[0]);
    out.output_root = nuxsec::app::trim(args[1]);
    if (args.size() == 3)
    {
        out.nthreads = std::stoi(nuxsec::app::trim(args[2]));
    }

    if (out.list_path.empty() || out.output_root.empty())
    {
        throw std::runtime_error("Invalid arguments (empty path)");
    }
    if (out.nthreads < 1)
    {
        throw std::runtime_error("Invalid thread count (must be >= 1)");
    }

    return out;
}

inline int run_analysis(const AnalysisArgs &analysis_args, const std::string &log_prefix)
{
    if (analysis_args.nthreads > 1)
    {
        ROOT::EnableImplicitMT(analysis_args.nthreads);
    }
    TH1::SetDefaultSumw2(true);

    const auto &analysis = nuxsec::AnalysisConfigService::instance();
    const auto entries = nuxsec::app::read_sample_list(analysis_args.list_path);
    const auto &specs = analysis.templates_1d();

    struct AnalysisInput
    {
        nuxsec::app::SampleListEntry entry;
        sample::SampleIO::Sample sample;
    };

    std::vector<nuxsec::io::AnalysisSampleRef> sample_refs;
    sample_refs.reserve(entries.size());
    std::vector<AnalysisInput> inputs;
    inputs.reserve(entries.size());

    for (const auto &entry : entries)
    {
        sample::SampleIO::Sample sample = sample::SampleIO::read(entry.output_path);

        nuxsec::io::AnalysisSampleRef ref;
        ref.sample_name = sample.sample_name;
        ref.sample_rootio_path = entry.output_path;
        ref.sample_kind = static_cast<int>(sample.kind);
        ref.beam_mode = static_cast<int>(sample.beam);
        ref.subrun_pot_sum = sample.subrun_pot_sum;
        ref.db_tortgt_pot_sum = sample.db_tortgt_pot_sum;
        ref.db_tor101_pot_sum = sample.db_tor101_pot_sum;
        sample_refs.push_back(ref);

        AnalysisInput input;
        input.entry = entry;
        input.sample = std::move(sample);
        inputs.push_back(std::move(input));
    }

    nuxsec::io::AnalysisHeader header;
    header.analysis_name = analysis.name();
    header.analysis_tree = analysis.tree_name();
    header.sample_list_source = analysis_args.list_path;

    nuxsec::io::AnalysisIO::init(analysis_args.output_root,
                                header,
                                sample_refs,
                                analysis.templates_1d_to_tsv(),
                                "compiled");

    nuxsec::io::AnalysisIO analysis_io(analysis_args.output_root, nuxsec::io::AnalysisIO::OpenMode::kUpdate);

    for (const auto &input : inputs)
    {
        const sample::SampleIO::Sample &sample = input.sample;
        ROOT::RDataFrame rdf = nuxsec::RDataFrameFactory::load_sample(sample, analysis.tree_name());
        const nuxsec::ProcessorEntry proc_entry = analysis.make_processor_entry(sample);

        const auto &processor = nuxsec::ColumnDerivationService::instance();
        ROOT::RDF::RNode node = processor.define(rdf, proc_entry);

        std::vector<ROOT::RDF::RResultPtr<TH1D>> booked;
        booked.reserve(specs.size());
        std::vector<ROOT::RDF::RResultHandle> handles;
        handles.reserve(specs.size());

        for (const auto &spec : specs)
        {
            ROOT::RDF::RNode filtered = node;
            if (!spec.selection.empty())
            {
                filtered = node.Filter(spec.selection, spec.name);
            }
            const std::string weight = spec.weight.empty() ? "w_template" : spec.weight;
            ROOT::RDF::TH1DModel model(spec.name.c_str(),
                                       spec.title.c_str(),
                                       spec.nbins,
                                       spec.xmin,
                                       spec.xmax);
            auto hist = filtered.Histo1D(model, spec.variable, weight);
            booked.push_back(hist);
            handles.emplace_back(hist);
        }

        ROOT::RDF::RunGraphs(handles);

        std::vector<std::pair<std::string, const TH1 *>> to_write;
        to_write.reserve(specs.size());
        for (size_t i = 0; i < specs.size(); ++i)
        {
            const TH1D &hist = booked[i].GetValue();
            to_write.emplace_back(specs[i].name, &hist);
        }

        analysis_io.put_histograms(nuxsec::io::AnalysisIO::kFamilyTemplates1D, sample.sample_name, to_write);

        std::cerr << "[" << log_prefix << "] analysis=" << analysis.name()
                  << " sample=" << sample.sample_name
                  << " kind=" << sample::SampleIO::sample_kind_name(sample.kind)
                  << " beam=" << sample::SampleIO::beam_mode_name(sample.beam)
                  << " templates=" << specs.size()
                  << " output=" << analysis_args.output_root
                  << "\n";
    }

    analysis_io.flush();

    return 0;
}

} // namespace app

} // namespace nuxsec

#endif // NUXSEC_APPS_ANALYSISCOMMAND_H
