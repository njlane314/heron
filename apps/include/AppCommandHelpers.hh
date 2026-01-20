/* -- C++ -- */
#ifndef Nuxsec_APPS_APPCOMMANDHELPERS_H_INCLUDED
#define Nuxsec_APPS_APPCOMMANDHELPERS_H_INCLUDED

#include <filesystem>
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
#include "ArtFileProvenanceRootIO.hh"
#include "RDataFrameFactory.hh"
#include "RunInfoSqliteReader.hh"
#include "SampleAggregator.hh"
#include "SampleRootIO.hh"
#include "SampleTypes.hh"
#include "SubrunTreeScanner.hh"
#include "TemplateRootIO.hh"

namespace nuxsec
{

namespace app
{

inline bool is_selection_data_file(const std::string &path)
{
    const auto pos = path.find_last_of("/\\");
    const std::string name = (pos == std::string::npos) ? path : path.substr(pos + 1);
    return name == "nuselection_data.root";
}

struct ArtArgs
{
    std::string artio_path;
    nuxsec::StageCfg stage_cfg;
    nuxsec::SampleKind sample_kind = nuxsec::SampleKind::kUnknown;
    nuxsec::BeamMode beam_mode = nuxsec::BeamMode::kUnknown;
};

inline ArtArgs parse_art_spec(const std::string &spec)
{
    std::vector<std::string> fields;
    size_t start = 0;
    while (start <= spec.size())
    {
        const size_t pos = spec.find(':', start);
        if (pos == std::string::npos)
        {
            fields.push_back(nuxsec::app::trim(spec.substr(start)));
            break;
        }
        fields.push_back(nuxsec::app::trim(spec.substr(start, pos - start)));
        start = pos + 1;
    }

    if (fields.size() < 2)
    {
        throw std::runtime_error("Bad stage spec (expected NAME:FILELIST): " + spec);
    }

    ArtArgs out;
    out.stage_cfg.stage_name = fields[0];
    out.stage_cfg.filelist_path = fields[1];

    if (out.stage_cfg.stage_name.empty() || out.stage_cfg.filelist_path.empty())
    {
        throw std::runtime_error("Bad stage spec: " + spec);
    }

    if (fields.size() >= 4)
    {
        out.sample_kind = nuxsec::parse_sample_kind(fields[2]);
        out.beam_mode = nuxsec::parse_beam_mode(fields[3]);
        if (out.sample_kind == nuxsec::SampleKind::kUnknown)
        {
            throw std::runtime_error("Bad stage sample kind: " + fields[2]);
        }
        if (out.beam_mode == nuxsec::BeamMode::kUnknown)
        {
            throw std::runtime_error("Bad stage beam mode: " + fields[3]);
        }
    }
    else if (fields.size() != 2)
    {
        throw std::runtime_error("Bad stage spec (expected NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]): " + spec);
    }

    out.artio_path = "build/out/art/art_prov_" + out.stage_cfg.stage_name + ".root";

    return out;
}

inline ArtArgs parse_art_args(const std::vector<std::string> &args, const std::string &usage)
{
    if (args.size() != 1)
    {
        throw std::runtime_error(usage);
    }

    return parse_art_spec(args[0]);
}

inline int run_artio(const ArtArgs &art_args, const std::string &log_prefix)
{
    const std::string db_path = "/exp/uboone/data/uboonebeam/beamdb/run.db";
    const double pot_scale = 1e12;

    std::filesystem::path out_path(art_args.artio_path);
    if (!out_path.parent_path().empty())
    {
        std::filesystem::create_directories(out_path.parent_path());
    }

    nuxsec::RunInfoSqliteReader db(db_path);

    const auto files = nuxsec::app::read_file_list(art_args.stage_cfg.filelist_path);

    nuxsec::ArtFileProvenance rec;
    rec.cfg = art_args.stage_cfg;
    rec.input_files = files;
    rec.kind = art_args.sample_kind;
    rec.beam = art_args.beam_mode;

    if (rec.kind == nuxsec::SampleKind::kUnknown && is_selection_data_file(files.front()))
    {
        rec.kind = nuxsec::SampleKind::kData;
    }

    rec.subrun = nuxsec::scanSubrunTree(files);
    rec.runinfo = db.sumRunInfo(rec.subrun.unique_pairs);

    rec.subrun.pot_sum *= pot_scale;
    rec.runinfo.tortgt_sum *= pot_scale;
    rec.runinfo.tor101_sum *= pot_scale;
    rec.runinfo.tor860_sum *= pot_scale;
    rec.runinfo.tor875_sum *= pot_scale;
    rec.db_tortgt_pot = rec.runinfo.tortgt_sum;
    rec.db_tor101_pot = rec.runinfo.tor101_sum;
    rec.scale = pot_scale;

    std::cerr << "[" << log_prefix << "] add stage=" << rec.cfg.stage_name
              << " files=" << rec.input_files.size()
              << " pairs=" << rec.subrun.unique_pairs.size()
              << " pot_sum=" << rec.subrun.pot_sum
              << " tortgt=" << rec.runinfo.tortgt_sum
              << "\n";

    nuxsec::ArtFileProvenanceRootIO::write(rec, art_args.artio_path);

    return 0;
}

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

    nuxsec::Sample sample = nuxsec::SampleAggregator::aggregate(sample_args.sample_name, files);
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

    nuxsec::TemplateRootIO::write_string_meta(tpl_args.output_root, "__global__", "analysis_name",
                                              analysis.Name());
    nuxsec::TemplateRootIO::write_string_meta(tpl_args.output_root,
                                              "__global__",
                                              "analysis_tree",
                                              analysis.TreeName());
    nuxsec::TemplateRootIO::write_string_meta(tpl_args.output_root,
                                              "__global__",
                                              "template_spec_source",
                                              "compiled");
    nuxsec::TemplateRootIO::write_string_meta(tpl_args.output_root,
                                              "__global__",
                                              "template_spec_tsv",
                                              analysis.Templates1DToTsv());

    for (const auto &entry : entries)
    {
        const nuxsec::Sample sample = nuxsec::SampleRootIO::read(entry.output_path);
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

        nuxsec::TemplateRootIO::write_histograms(tpl_args.output_root, sample.sample_name, to_write);
        nuxsec::TemplateRootIO::write_string_meta(tpl_args.output_root,
                                                  sample.sample_name,
                                                  "sample_kind",
                                                  nuxsec::sample_kind_name(sample.kind));
        nuxsec::TemplateRootIO::write_string_meta(tpl_args.output_root,
                                                  sample.sample_name,
                                                  "beam_mode",
                                                  nuxsec::beam_mode_name(sample.beam));
        nuxsec::TemplateRootIO::write_string_meta(tpl_args.output_root,
                                                  sample.sample_name,
                                                  "sample_rootio_path",
                                                  entry.output_path);

        std::cerr << "[" << log_prefix << "] analysis=" << analysis.Name()
                  << " sample=" << sample.sample_name
                  << " kind=" << nuxsec::sample_kind_name(sample.kind)
                  << " beam=" << nuxsec::beam_mode_name(sample.beam)
                  << " templates=" << specs.size()
                  << " output=" << tpl_args.output_root
                  << "\n";
    }

    return 0;
}

}

}

#endif
