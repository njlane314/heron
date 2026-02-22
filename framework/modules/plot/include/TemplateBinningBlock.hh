/* -- C++ -- */
/**
 *  @file framework/modules/plot/include/TemplateBinningBlock.hh
 *
 *  @brief Block-style binning container for template-based optimisation output.
 */

#ifndef heron_plot_template_binning_block_H
#define heron_plot_template_binning_block_H

#include <string>
#include <vector>

#include <TNamed.h>

class TemplateBinningBlock : public TNamed
{
  public:
    TemplateBinningBlock(const std::string &name,
                         const std::string &title,
                         const std::vector<double> &edges,
                         const std::string &selection,
                         int bin_type);
    virtual ~TemplateBinningBlock() = default;

    int GetNBinsX() const;
    int GetNBinsY(int idx) const;
    int GetNBinsZ(int idx, int idy) const;

    std::string GetBinDef(int idx) const;
    std::string GetBinDef(int x, int y) const;
    std::string GetBinDef(int x, int y, int z) const;

    int GetBinType() const;

    Double_t GetBinXLow(int i) const;
    Double_t GetBinXHigh(int i) const;
    Double_t GetBinYLow(int i, int j) const;
    Double_t GetBinYHigh(int i, int j) const;

    std::string GetXName() const;
    std::string GetXNameUnit() const;
    std::string GetXTitle() const;
    std::string GetXTitleUnit() const;
    std::string GetXTexTitle() const;
    std::string GetXTexTitleUnit() const;

    std::string GetYName() const;
    std::string GetYNameUnit() const;
    std::string GetYTitle() const;
    std::string GetYTitleUnit() const;
    std::string GetYTexTitle() const;
    std::string GetYTexTitleUnit() const;

    std::string GetZName() const;
    std::string GetZNameUnit() const;
    std::string GetZTitle() const;
    std::string GetZTitleUnit() const;
    std::string GetZTexTitle() const;
    std::string GetZTexTitleUnit() const;

    std::string GetSelection() const;

    std::vector<double> GetVector() const;
    std::vector<double> GetVector(int x) const;
    std::vector<double> GetVector(int x, int y) const;

    bool Is1D() const;
    bool Is2D() const;
    bool Is3D() const;

  private:
    void init();

    std::string m_name;
    std::string m_title;
    std::vector<double> m_edges;
    std::string m_selection;
    int m_bin_type = 0;

    std::vector<std::string> m_bin_def;

    std::string m_x_name;
    std::string m_x_name_unit;
    std::string m_x_title;
    std::string m_x_title_unit;
    std::string m_x_tex_title;
    std::string m_x_tex_title_unit;

    std::string m_y_name;
    std::string m_y_name_unit;
    std::string m_y_title;
    std::string m_y_title_unit;
    std::string m_y_tex_title;
    std::string m_y_tex_title_unit;

    std::string m_z_name;
    std::string m_z_name_unit;
    std::string m_z_title;
    std::string m_z_title_unit;
    std::string m_z_tex_title;
    std::string m_z_tex_title_unit;
};

class TemplateBinningBlockReco
{
  public:
    explicit TemplateBinningBlockReco(TemplateBinningBlock *p_block_reco)
        : p_block_reco_(p_block_reco) {}
    virtual ~TemplateBinningBlockReco() = default;

    TemplateBinningBlock *p_block_reco_ = nullptr;
};

#endif // heron_plot_template_binning_block_H
