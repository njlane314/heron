/* -- C++ -- */
#ifndef Nuxsec_APPS_APPTEMPLATECOMMAND_H_INCLUDED
#define Nuxsec_APPS_APPTEMPLATECOMMAND_H_INCLUDED

#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TH1.h>
#include <TH1D.h>

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>

#include "AnalysisDefinition.hh"
#include "AnalysisRdfDefinitions.hh"
#include "AppUtils.hh"
#include "RDataFrameFactory.hh"
#include "SampleIO.hh"
#include "TemplateIO.hh"

namespace nuxsec
{

namespace app
{

struct TemplateArgs
{
    std::string list_path;
    std::string output_root;
    int nthreads = 1;
};

inline TemplateArgs parse_template_args(const std::vector<std::string> &args, const std::string &usage)
{
    if (args.size() < 2 || args.size() > 3)
    {
        throw std::runtime_error(usage);
    }

    TemplateArgs out;
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

inline int run_template(const TemplateArgs &tpl_args, const std::string &log_prefix)
{
    if (tpl_args.nthreads > 1)
    {
        ROOT::EnableImplicitMT(tpl_args.nthreads);
    }
    TH1::SetDefaultSumw2(true);

    const auto &analysis = nuxsec::AnalysisDefinition::Instance();
    const auto entries = nuxsec::app::read_sample_list(tpl_args.list_path);
    const auto &specs = analysis.Templates1D();

    nuxsec::TemplateIO::write_string_meta(tpl_args.output_root, "__global__", "analysis_name",
                                              analysis.Name());
    nuxsec::TemplateIO::write_string_meta(tpl_args.output_root,
                                              "__global__",
                                              "analysis_tree",
                                              analysis.TreeName());
    nuxsec::TemplateIO::write_string_meta(tpl_args.output_root,
                                              "__global__",
                                              "template_spec_source",
                                              "compiled");
    nuxsec::TemplateIO::write_string_meta(tpl_args.output_root,
                                              "__global__",
                                              "template_spec_tsv",
                                              analysis.Templates1DToTsv());

    for (const auto &entry : entries)
    {
        const sample::SampleIO::Sample sample = sample::SampleIO::read(entry.output_path);
        ROOT::RDataFrame rdf = nuxsec::RDataFrameFactory::load_sample(sample, analysis.TreeName());
        const nuxsec::ProcessorEntry proc_entry = analysis.MakeProcessorEntry(sample);

        const auto &processor = nuxsec::AnalysisRdfDefinitions::Instance();
        ROOT::RDF::RNode node = processor.Define(rdf, proc_entry);

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

        nuxsec::TemplateIO::write_histograms(tpl_args.output_root, sample.sample_name, to_write);
        nuxsec::TemplateIO::write_string_meta(tpl_args.output_root,
                                                  sample.sample_name,
                                                  "sample_kind",
                                                  sample::SampleIO::sample_kind_name(sample.kind));
        nuxsec::TemplateIO::write_string_meta(tpl_args.output_root,
                                                  sample.sample_name,
                                                  "beam_mode",
                                                  sample::SampleIO::beam_mode_name(sample.beam));
        nuxsec::TemplateIO::write_string_meta(tpl_args.output_root,
                                                  sample.sample_name,
                                                  "sample_rootio_path",
                                                  entry.output_path);

        std::cerr << "[" << log_prefix << "] analysis=" << analysis.Name()
                  << " sample=" << sample.sample_name
                  << " kind=" << sample::SampleIO::sample_kind_name(sample.kind)
                  << " beam=" << sample::SampleIO::beam_mode_name(sample.beam)
                  << " templates=" << specs.size()
                  << " output=" << tpl_args.output_root
                  << "\n";
    }

    return 0;
}

}

}

#endif
