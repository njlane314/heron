/* -- C++ -- */
/**
 *  @file  io/src/TemplateRootIO.cc
 *
 *  @brief Implementation for TemplateRootIO helpers.
 */

#include "TemplateRootIO.hh"

#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TNamed.h>
#include <TParameter.h>

#include <memory>
#include <stdexcept>

namespace
{

TDirectory *get_or_make_dir(TDirectory *parent, const std::string &name)
{
    if (!parent)
    {
        return nullptr;
    }
    if (name.empty())
    {
        return parent;
    }
    if (auto *existing = parent->GetDirectory(name.c_str()))
    {
        return existing;
    }
    return parent->mkdir(name.c_str());
}

}

namespace nuxsec
{

void TemplateRootIO::write_histograms(const std::string &root_path,
                                      const std::string &sample_name,
                                      const std::vector<std::pair<std::string, const TH1 *>> &hists,
                                      const TemplateWriteOptions &opt)
{
    if (hists.empty())
    {
        return;
    }

    std::unique_ptr<TFile> fout(TFile::Open(root_path.c_str(), "UPDATE"));
    if (!fout || fout->IsZombie())
    {
        throw std::runtime_error("TemplateRootIO::write_histograms: cannot open " + root_path);
    }

    TDirectory *top = get_or_make_dir(fout.get(), opt.top_dir);
    TDirectory *sdir = get_or_make_dir(top, sample_name);
    TDirectory *hdir = get_or_make_dir(sdir, "hists");
    if (!hdir)
    {
        throw std::runtime_error("TemplateRootIO::write_histograms: failed to create directories");
    }

    hdir->cd();
    for (const auto &entry : hists)
    {
        const std::string &hname = entry.first;
        const TH1 *hptr = entry.second;
        if (!hptr)
        {
            continue;
        }

        std::unique_ptr<TH1> clone(dynamic_cast<TH1 *>(hptr->Clone(hname.c_str())));
        if (!clone)
        {
            throw std::runtime_error("TemplateRootIO::write_histograms: failed to clone TH1 " + hname);
        }

        clone->SetDirectory(hdir);
        const int wopt = opt.overwrite ? TObject::kOverwrite : 0;
        clone->Write(hname.c_str(), wopt);
    }

    fout->Write();
    fout->Close();
}

void TemplateRootIO::write_string_meta(const std::string &root_path,
                                       const std::string &sample_name,
                                       const std::string &key,
                                       const std::string &value,
                                       const TemplateWriteOptions &opt)
{
    std::unique_ptr<TFile> fout(TFile::Open(root_path.c_str(), "UPDATE"));
    if (!fout || fout->IsZombie())
    {
        throw std::runtime_error("TemplateRootIO::write_string_meta: cannot open " + root_path);
    }

    TDirectory *top = get_or_make_dir(fout.get(), opt.top_dir);
    TDirectory *sdir = get_or_make_dir(top, sample_name);
    TDirectory *mdir = get_or_make_dir(sdir, "meta");
    if (!mdir)
    {
        throw std::runtime_error("TemplateRootIO::write_string_meta: failed to create meta directory");
    }

    mdir->cd();
    TNamed named(key.c_str(), value.c_str());
    const int wopt = opt.overwrite ? TObject::kOverwrite : 0;
    named.Write(key.c_str(), wopt);

    fout->Write();
    fout->Close();
}

void TemplateRootIO::write_double_meta(const std::string &root_path,
                                       const std::string &sample_name,
                                       const std::string &key,
                                       double value,
                                       const TemplateWriteOptions &opt)
{
    std::unique_ptr<TFile> fout(TFile::Open(root_path.c_str(), "UPDATE"));
    if (!fout || fout->IsZombie())
    {
        throw std::runtime_error("TemplateRootIO::write_double_meta: cannot open " + root_path);
    }

    TDirectory *top = get_or_make_dir(fout.get(), opt.top_dir);
    TDirectory *sdir = get_or_make_dir(top, sample_name);
    TDirectory *mdir = get_or_make_dir(sdir, "meta");
    if (!mdir)
    {
        throw std::runtime_error("TemplateRootIO::write_double_meta: failed to create meta directory");
    }

    mdir->cd();
    TParameter<double> par(key.c_str(), value);
    const int wopt = opt.overwrite ? TObject::kOverwrite : 0;
    par.Write(key.c_str(), wopt);

    fout->Write();
    fout->Close();
}

void TemplateRootIO::write_syst_histograms(const std::string &root_path,
                                           const std::string &sample_name,
                                           const std::string &syst_name,
                                           const std::string &variation,
                                           const std::vector<std::pair<std::string, const TH1 *>> &hists,
                                           const TemplateWriteOptions &opt)
{
    if (hists.empty())
    {
        return;
    }
    if (syst_name.empty() || variation.empty())
    {
        return;
    }

    std::unique_ptr<TFile> fout(TFile::Open(root_path.c_str(), "UPDATE"));
    if (!fout || fout->IsZombie())
    {
        throw std::runtime_error("TemplateRootIO::write_syst_histograms: cannot open " + root_path);
    }

    TDirectory *top = get_or_make_dir(fout.get(), opt.top_dir);
    TDirectory *sdir = get_or_make_dir(top, sample_name);
    TDirectory *sysd = get_or_make_dir(sdir, "systs");
    TDirectory *thiss = get_or_make_dir(sysd, syst_name);
    TDirectory *vard = get_or_make_dir(thiss, variation);
    TDirectory *hdir = get_or_make_dir(vard, "hists");
    if (!hdir)
    {
        throw std::runtime_error("TemplateRootIO::write_syst_histograms: failed to create directories");
    }

    hdir->cd();
    for (const auto &entry : hists)
    {
        const std::string &hname = entry.first;
        const TH1 *hptr = entry.second;
        if (!hptr)
        {
            continue;
        }

        std::unique_ptr<TH1> clone(dynamic_cast<TH1 *>(hptr->Clone(hname.c_str())));
        if (!clone)
        {
            throw std::runtime_error("TemplateRootIO::write_syst_histograms: failed to clone TH1 " + hname);
        }

        clone->SetDirectory(hdir);
        const int wopt = opt.overwrite ? TObject::kOverwrite : 0;
        clone->Write(hname.c_str(), wopt);
    }

    fout->Write();
    fout->Close();
}

void TemplateRootIO::write_syst_flag_meta(const std::string &root_path,
                                          const std::string &syst_name,
                                          const std::string &key,
                                          const std::string &value,
                                          const TemplateWriteOptions &opt)
{
    if (syst_name.empty())
    {
        return;
    }

    std::unique_ptr<TFile> fout(TFile::Open(root_path.c_str(), "UPDATE"));
    if (!fout || fout->IsZombie())
    {
        throw std::runtime_error("TemplateRootIO::write_syst_flag_meta: cannot open " + root_path);
    }

    TDirectory *top = get_or_make_dir(fout.get(), opt.top_dir);
    TDirectory *gdir = get_or_make_dir(top, "__global__");
    TDirectory *mdir = get_or_make_dir(gdir, "meta");
    TDirectory *sdir = get_or_make_dir(mdir, "systs");
    TDirectory *thiss = get_or_make_dir(sdir, syst_name);
    if (!thiss)
    {
        throw std::runtime_error("TemplateRootIO::write_syst_flag_meta: failed to create syst meta directories");
    }

    thiss->cd();
    TNamed named(key.c_str(), value.c_str());
    const int wopt = opt.overwrite ? TObject::kOverwrite : 0;
    named.Write(key.c_str(), wopt);

    fout->Write();
    fout->Close();
}

} // namespace nuxsec
