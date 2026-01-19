/* -- C++ -- */
/**
 *  @file  apps/src/nuxsec.cc
 *
 *  @brief Unified CLI for Nuxsec utilities.
 */

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <TROOT.h>

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

namespace
{

const char *kUsageArt = "Usage: nuxsec a|art|artio|artio-aggregate NAME:FILELIST[:SAMPLE_KIND:BEAM_MODE]";
const char *kUsageSample = "Usage: nuxsec s|samp|sample|sample-aggregate NAME:FILELIST";
const char *kUsageTemplate = "Usage: nuxsec t|tpl|template|template-make SAMPLE_LIST.tsv OUTPUT.root [NTHREADS]";
const char *kUsagePlot = "Usage: nuxsec plot [OUTSTEM]";

bool is_help_arg(const std::string &arg)
{
    return arg == "-h" || arg == "--help";
}

void print_main_help(std::ostream &out)
{
    out << "Usage: nuxsec <command> [args]\n\n"
        << "Commands:\n"
        << "  a|art|artio     (artio-aggregate)  Aggregate art provenance for a stage\n"
        << "  s|samp|sample   (sample-aggregate) Aggregate Sample ROOT files from art provenance\n"
        << "  t|tpl|template  (template-make)    Build template histograms from sample list\n"
        << "  plot                         Build POT timeline plot via ROOT macro\n"
        << "\nRun 'nuxsec <command> --help' for command-specific usage.\n";
}

bool is_selection_data_file(const std::string &path)
{
    const auto pos = path.find_last_of("/\\");
    const std::string name = (pos == std::string::npos) ? path : path.substr(pos + 1);
    return name == "nuselection_data.root";
}

std::filesystem::path find_repo_root()
{
    std::vector<std::filesystem::path> candidates;

    std::error_code ec;
    const auto exe = std::filesystem::read_symlink("/proc/self/exe", ec);
    if (!ec)
    {
        candidates.push_back(exe.parent_path());
    }
    candidates.push_back(std::filesystem::current_path());

    for (auto base : candidates)
    {
        for (int i = 0; i < 6; ++i)
        {
            if (std::filesystem::exists(base / "plot/macro/NuxsecPlotDriver.C"))
            {
                return base;
            }
            if (!base.has_parent_path())
            {
                break;
            }
            base = base.parent_path();
        }
    }

    return std::filesystem::current_path();
}


struct ArtArgs
{
    std::string artio_path;
    nuxsec::StageCfg stage_cfg;
    nuxsec::SampleKind sample_kind = nuxsec::SampleKind::kUnknown;
    nuxsec::BeamMode beam_mode = nuxsec::BeamMode::kUnknown;
};

ArtArgs parse_art_args(const std::vector<std::string> &args)
{
    if (args.size() != 1)
    {
        throw std::runtime_error(kUsageArt);
    }

    const std::string spec = args[0];

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

int run_artio(const std::vector<std::string> &args)
{
    if (args.empty() || (args.size() == 1 && is_help_arg(args[0])))
    {
        std::cout << kUsageArt << "\n";
        return 0;
    }

    const std::string db_path = "/exp/uboone/data/uboonebeam/beamdb/run.db";
    const double pot_scale = 1e12;

    const ArtArgs art_args = parse_art_args(args);
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

    std::cerr << "[nuxsec artio-aggregate] add stage=" << rec.cfg.stage_name
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

SampleArgs parse_sample_args(const std::vector<std::string> &args)
{
    if (args.size() != 1)
    {
        throw std::runtime_error(kUsageSample);
    }

    const std::string spec = args[0];
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
        if (errno == ENOENT)
        {
            return {};
        }
        throw std::runtime_error("Failed to open sample list: " + list_path +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }

    std::vector<SampleListEntry> entries;
    std::string line;
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
        SampleListEntry entry;
        entry.sample_name = fields[0];
        entry.sample_kind = fields[1];
        entry.beam_mode = fields[2];
        entry.output_path = fields[3];
        entries.push_back(std::move(entry));
    }
    return entries;
}

void write_sample_list(const std::string &list_path, std::vector<SampleListEntry> entries)
{
    std::sort(entries.begin(), entries.end(),
              [](const SampleListEntry &a, const SampleListEntry &b)
              {
                  return std::tie(a.sample_kind, a.beam_mode, a.sample_name) <
                         std::tie(b.sample_kind, b.beam_mode, b.sample_name);
              });

    std::ofstream fout(list_path, std::ios::trunc);
    if (!fout)
    {
        throw std::runtime_error("Failed to open sample list for writing: " + list_path +
                                 " (errno=" + std::to_string(errno) + " " + std::strerror(errno) + ")");
    }
    fout << "# sample_name\tsample_kind\tbeam_mode\toutput_path\n";
    for (const auto &entry : entries)
    {
        fout << entry.sample_name << "\t"
             << entry.sample_kind << "\t"
             << entry.beam_mode << "\t"
             << entry.output_path << "\n";
    }
}

void update_sample_list(const std::string &list_path, const nuxsec::Sample &sample,
                        const std::string &output_path)
{
    auto entries = read_sample_list(list_path);
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
        SampleListEntry entry;
        entry.sample_name = sample.sample_name;
        entry.sample_kind = kind_name;
        entry.beam_mode = beam_name;
        entry.output_path = output_path;
        entries.push_back(std::move(entry));
    }

    write_sample_list(list_path, std::move(entries));
}

int run_sample(const std::vector<std::string> &args)
{
    if (args.empty() || (args.size() == 1 && is_help_arg(args[0])))
    {
        std::cout << kUsageSample << "\n";
        return 0;
    }

    const SampleArgs sample_args = parse_sample_args(args);
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

    std::cerr << "[nuxsec sample-aggregate] sample=" << sample.sample_name
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

struct TemplateListEntry
{
    std::string sample_name;
    std::string sample_kind;
    std::string beam_mode;
    std::string output_path;
};

std::vector<TemplateListEntry> read_template_list(const std::string &list_path)
{
    std::ifstream fin(list_path);
    if (!fin)
    {
        throw std::runtime_error("Failed to open sample list: " + list_path);
    }

    std::vector<TemplateListEntry> entries;
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

        TemplateListEntry entry;
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

struct TemplateArgs
{
    std::string list_path;
    std::string output_root;
    int nthreads = 1;
};

TemplateArgs parse_template_args(const std::vector<std::string> &args)
{
    if (args.size() < 2 || args.size() > 3)
    {
        throw std::runtime_error(kUsageTemplate);
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

int run_template(const std::vector<std::string> &args)
{
    if (args.empty() || (args.size() == 1 && is_help_arg(args[0])))
    {
        std::cout << kUsageTemplate << "\n";
        return 0;
    }

    const TemplateArgs tpl_args = parse_template_args(args);
    if (tpl_args.nthreads > 1)
    {
        ROOT::EnableImplicitMT(tpl_args.nthreads);
    }
    TH1::SetDefaultSumw2(true);

    const auto &analysis = nuxsec::AnalysisDefinition::Instance();
    const auto entries = read_template_list(tpl_args.list_path);
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

        std::cerr << "[nuxsec template-make] analysis=" << analysis.Name()
                  << " sample=" << sample.sample_name
                  << " kind=" << nuxsec::sample_kind_name(sample.kind)
                  << " beam=" << nuxsec::beam_mode_name(sample.beam)
                  << " templates=" << specs.size()
                  << " output=" << tpl_args.output_root
                  << "\n";
    }

    return 0;
}

struct PlotArgs
{
    std::string outstem = "build/out/plot/pot_timeline";
};

PlotArgs parse_plot_args(const std::vector<std::string> &args)
{
    if (args.size() > 1)
    {
        throw std::runtime_error(kUsagePlot);
    }

    PlotArgs out;
    if (args.size() == 1)
    {
        out.outstem = nuxsec::app::trim(args[0]);
    }

    if (out.outstem.empty())
    {
        throw std::runtime_error("Invalid output stem");
    }
    return out;
}

int run_plot(const std::vector<std::string> &args)
{
    if (args.empty() || (args.size() == 1 && is_help_arg(args[0])))
    {
        std::cout << kUsagePlot << "\n";
        return 0;
    }

    const PlotArgs plot_args = parse_plot_args(args);
    std::filesystem::path outstem_path(plot_args.outstem);
    if (!outstem_path.parent_path().empty())
    {
        std::filesystem::create_directories(outstem_path.parent_path());
    }
    const auto repo_root = find_repo_root();
    const auto macro_path = repo_root / "plot/macro/NuxsecPlotDriver.C";
    if (!std::filesystem::exists(macro_path))
    {
        throw std::runtime_error("Plot macro not found at " + macro_path.string());
    }

    const std::string load_cmd = ".L " + macro_path.string();
    gROOT->ProcessLine(load_cmd.c_str());

    const std::string call_cmd = "nuxsec_plot(\"" + plot_args.outstem + "\")";
    const long result = gROOT->ProcessLine(call_cmd.c_str());
    return static_cast<int>(result);
}

}

int main(int argc, char **argv)
{
    try
    {
        if (argc < 2)
        {
            print_main_help(std::cerr);
            return 1;
        }

        const std::string command = argv[1];
        std::vector<std::string> args;
        args.reserve(static_cast<size_t>(argc - 2));
        for (int i = 2; i < argc; ++i)
        {
            args.emplace_back(argv[i]);
        }

        if (command == "help" || command == "-h" || command == "--help")
        {
            print_main_help(std::cout);
            return 0;
        }

        if (command == "artio-aggregate" || command == "artio" || command == "art" || command == "a")
        {
            return run_artio(args);
        }
        if (command == "sample-aggregate" || command == "sample" || command == "samp" || command == "s")
        {
            return run_sample(args);
        }
        if (command == "template-make" || command == "template" || command == "tpl" || command == "t")
        {
            return run_template(args);
        }
        if (command == "plot")
        {
            return run_plot(args);
        }

        std::cerr << "Unknown command: " << command << "\n";
        print_main_help(std::cerr);
        return 1;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
