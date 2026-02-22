/* -- C++ -- */
/**
 *  @file framework/modules/plot/src/TemplateBinningBlock.cc
 *
 *  @brief Implementation of block-style 1D binning container.
 */

#include "TemplateBinningBlock.hh"

#include <sstream>
#include <stdexcept>

TemplateBinningBlock::TemplateBinningBlock(const std::string &name,
                                           const std::string &title,
                                           const std::vector<double> &edges,
                                           const std::string &selection,
                                           int bin_type)
    : TNamed(name.c_str(), title.c_str()),
      m_name(name),
      m_title(title),
      m_edges(edges),
      m_selection(selection),
      m_bin_type(bin_type)
{
    init();
}

void TemplateBinningBlock::init()
{
    if (m_edges.size() < 2)
        throw std::runtime_error("TemplateBinningBlock requires at least two edges");

    m_bin_def.clear();
    for (size_t i = 0; i + 1 < m_edges.size(); ++i)
    {
        std::ostringstream os;
        os << "[" << m_edges[i] << ", " << m_edges[i + 1] << ")";
        m_bin_def.push_back(os.str());
    }

    m_x_name = m_name;
    m_x_title = m_title;
    m_x_tex_title = m_title;
}

int TemplateBinningBlock::GetNBinsX() const { return static_cast<int>(m_edges.size()) - 1; }
int TemplateBinningBlock::GetNBinsY(int) const { return 0; }
int TemplateBinningBlock::GetNBinsZ(int, int) const { return 0; }

std::string TemplateBinningBlock::GetBinDef(int idx) const { return m_bin_def.at(static_cast<size_t>(idx)); }
std::string TemplateBinningBlock::GetBinDef(int x, int) const { return GetBinDef(x); }
std::string TemplateBinningBlock::GetBinDef(int x, int, int) const { return GetBinDef(x); }

int TemplateBinningBlock::GetBinType() const { return m_bin_type; }

Double_t TemplateBinningBlock::GetBinXLow(int i) const { return m_edges.at(static_cast<size_t>(i)); }
Double_t TemplateBinningBlock::GetBinXHigh(int i) const { return m_edges.at(static_cast<size_t>(i + 1)); }
Double_t TemplateBinningBlock::GetBinYLow(int, int) const { return -9999999; }
Double_t TemplateBinningBlock::GetBinYHigh(int, int) const { return -9999999; }

std::string TemplateBinningBlock::GetXName() const { return m_x_name; }
std::string TemplateBinningBlock::GetXNameUnit() const { return m_x_name_unit; }
std::string TemplateBinningBlock::GetXTitle() const { return m_x_title; }
std::string TemplateBinningBlock::GetXTitleUnit() const { return m_x_title_unit; }
std::string TemplateBinningBlock::GetXTexTitle() const { return m_x_tex_title; }
std::string TemplateBinningBlock::GetXTexTitleUnit() const { return m_x_tex_title_unit; }

std::string TemplateBinningBlock::GetYName() const { return m_y_name; }
std::string TemplateBinningBlock::GetYNameUnit() const { return m_y_name_unit; }
std::string TemplateBinningBlock::GetYTitle() const { return m_y_title; }
std::string TemplateBinningBlock::GetYTitleUnit() const { return m_y_title_unit; }
std::string TemplateBinningBlock::GetYTexTitle() const { return m_y_tex_title; }
std::string TemplateBinningBlock::GetYTexTitleUnit() const { return m_y_tex_title_unit; }

std::string TemplateBinningBlock::GetZName() const { return m_z_name; }
std::string TemplateBinningBlock::GetZNameUnit() const { return m_z_name_unit; }
std::string TemplateBinningBlock::GetZTitle() const { return m_z_title; }
std::string TemplateBinningBlock::GetZTitleUnit() const { return m_z_title_unit; }
std::string TemplateBinningBlock::GetZTexTitle() const { return m_z_tex_title; }
std::string TemplateBinningBlock::GetZTexTitleUnit() const { return m_z_tex_title_unit; }

std::string TemplateBinningBlock::GetSelection() const { return m_selection; }

std::vector<double> TemplateBinningBlock::GetVector() const { return m_edges; }
std::vector<double> TemplateBinningBlock::GetVector(int) const { return m_edges; }
std::vector<double> TemplateBinningBlock::GetVector(int, int) const { return m_edges; }

bool TemplateBinningBlock::Is1D() const { return true; }
bool TemplateBinningBlock::Is2D() const { return false; }
bool TemplateBinningBlock::Is3D() const { return false; }
