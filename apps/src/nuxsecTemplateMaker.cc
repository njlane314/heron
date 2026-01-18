/* -- C++ -- */
/**
 *  @file  apps/src/nuxsecTemplateMaker.cc
 *
 *  @brief Build binned template histograms from aggregated samples.
 */

#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TH1.h>
#include <TH1D.h>

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>

#include "AnalysisRdfDefinitions.hh"
#include "AppUtils.hh"
#include "RDataFrameFactory.hh"
#include "SampleRootIO.hh"
#include "SampleTypes.hh"
#include "TemplateRootIO.hh"
#include "TemplateSpec.hh"

namespace
{

struct SampleListEntry
{
    std::string sample_name;
    std::string sample_kind;
    std::string beam_mode;
    std::string output_path;
};

std::vector<SampleListEntry> read_sample_list(const std::string &list_path)
{
    std::ifstream fin(list_path);
    if (!fin)
    {
        throw std::runtime_error("Failed to open sample list: " + list_path);
    }

    std::vector<SampleListEntry> entries;
    std::string line;
    bool first_nonempty = true;
    while (std::getline(fin, line))
    {
        line = nuxsec::app::trim(line);
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        const auto fields = nuxsec::app::split_tabs(line);
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
        entry.sample_kind = fields[1];
        entry.beam_mode = fields[2];
        entry.output_path = fields[3];
        entries.push_back(std::move(entry));
    }

    if (entries.empty())
    {
        throw std::runtime_error("Sample list is empty: " + list_path);
    }

    return entries;
}

struct Args
{
    std::string list_path;
    std::string tree_name;
    std::string template_path;
    std::string output_root;
    int nthreads = 1;
};

Args parse_args(int argc, char **argv)
{
    if (argc < 5 || argc > 6)
    {
        throw std::runtime_error("Usage: nuxsecTemplateMaker SAMPLE_LIST.tsv TREE_NAME "
                                 "TEMPLATES.tsv OUTPUT.root [NTHREADS]");
    }

    Args args;
    args.list_path = nuxsec::app::trim(argv[1]);
    args.tree_name = nuxsec::app::trim(argv[2]);
    args.template_path = nuxsec::app::trim(argv[3]);
    args.output_root = nuxsec::app::trim(argv[4]);
    if (argc == 6)
    {
        args.nthreads = std::stoi(nuxsec::app::trim(argv[5]));
    }

    if (args.list_path.empty() || args.tree_name.empty() || args.template_path.empty() ||
        args.output_root.empty())
    {
        throw std::runtime_error("Invalid arguments (empty path or tree name)");
    }

    if (args.nthreads < 1)
    {
        throw std::runtime_error("Invalid thread count (must be >= 1)");
    }

    return args;
}

}

int main(int argc, char **argv)
{
    try
    {
        const Args args = parse_args(argc, argv);
        if (args.nthreads > 1)
        {
            ROOT::EnableImplicitMT(args.nthreads);
        }
        TH1::SetDefaultSumw2(true);

        const auto entries = read_sample_list(args.list_path);
        const auto specs = nuxsec::read_template_spec_1d_tsv(args.template_path);

        nuxsec::TemplateRootIO::write_string_meta(args.output_root,
                                                  "__global__",
                                                  "template_spec_path",
                                                  args.template_path);

        for (const auto &entry : entries)
        {
            const nuxsec::Sample sample = nuxsec::SampleRootIO::read(entry.output_path);
            ROOT::RDataFrame rdf = nuxsec::RDataFrameFactory::load_sample(sample, args.tree_name);

            nuxsec::ProcessorEntry proc_entry;
            switch (sample.kind)
            {
            case nuxsec::SampleKind::kData:
                proc_entry.source = nuxsec::SourceKind::kData;
                break;
            case nuxsec::SampleKind::kEXT:
                proc_entry.source = nuxsec::SourceKind::kExt;
                proc_entry.trig_nom = sample.db_tor101_pot_sum;
                proc_entry.trig_eqv = sample.subrun_pot_sum;
                break;
            case nuxsec::SampleKind::kOverlay:
            case nuxsec::SampleKind::kDirt:
            case nuxsec::SampleKind::kStrangeness:
                proc_entry.source = nuxsec::SourceKind::kMC;
                proc_entry.pot_nom = sample.db_tortgt_pot_sum;
                proc_entry.pot_eqv = sample.subrun_pot_sum;
                break;
            default:
                proc_entry.source = nuxsec::SourceKind::kUnknown;
                break;
            }

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

            nuxsec::TemplateRootIO::write_histograms(args.output_root, sample.sample_name, to_write);
            nuxsec::TemplateRootIO::write_string_meta(args.output_root,
                                                      sample.sample_name,
                                                      "sample_kind",
                                                      nuxsec::sample_kind_name(sample.kind));
            nuxsec::TemplateRootIO::write_string_meta(args.output_root,
                                                      sample.sample_name,
                                                      "beam_mode",
                                                      nuxsec::beam_mode_name(sample.beam));
            nuxsec::TemplateRootIO::write_string_meta(args.output_root,
                                                      sample.sample_name,
                                                      "sample_rootio_path",
                                                      entry.output_path);

            std::cerr << "[nuxsecTemplateMaker] sample=" << sample.sample_name
                      << " kind=" << entry.sample_kind
                      << " beam=" << entry.beam_mode
                      << " templates=" << specs.size()
                      << " output=" << args.output_root
                      << "\n";
        }

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
