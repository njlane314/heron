/* -- C++ -- */
/**
 *  @file  ana/src/TemplateSpec.cc
 *
 *  @brief Implementation of template specification parsing.
 */

#include "TemplateSpec.hh"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{

std::vector<std::string> split_tabs(const std::string &line)
{
    std::vector<std::string> out;
    std::string cell;
    std::stringstream ss(line);
    while (std::getline(ss, cell, '\t'))
    {
        out.push_back(cell);
    }
    return out;
}

std::string trim(std::string s)
{
    auto notspace = [](unsigned char c)
    {
        return std::isspace(c) == 0;
    };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notspace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
    return s;
}

bool looks_like_header(const std::vector<std::string> &cols)
{
    if (cols.empty())
    {
        return false;
    }
    return (cols[0] == "name" || cols[0] == "template_name");
}

}

namespace nuxsec
{

std::vector<TemplateSpec1D> read_template_spec_1d_tsv(const std::string &path)
{
    std::ifstream fin(path);
    if (!fin)
    {
        throw std::runtime_error("read_template_spec_1d_tsv: failed to open " + path);
    }

    std::vector<TemplateSpec1D> specs;
    std::string line;
    bool first_nonempty = true;
    while (std::getline(fin, line))
    {
        line = trim(line);
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        const auto cols = split_tabs(line);
        if (cols.size() < 8)
        {
            throw std::runtime_error("read_template_spec_1d_tsv: expected >=8 TSV columns, got " +
                                     std::to_string(cols.size()) + " in line: " + line);
        }

        if (first_nonempty && looks_like_header(cols))
        {
            first_nonempty = false;
            continue;
        }
        first_nonempty = false;

        TemplateSpec1D spec;
        spec.name = cols[0];
        spec.title = cols[1];
        spec.selection = cols[2];
        spec.variable = cols[3];
        spec.weight = cols[4];

        try
        {
            spec.nbins = std::stoi(cols[5]);
            spec.xmin = std::stod(cols[6]);
            spec.xmax = std::stod(cols[7]);
        }
        catch (const std::exception &)
        {
            throw std::runtime_error("read_template_spec_1d_tsv: bad numeric fields in line: " + line);
        }

        if (spec.title.empty())
        {
            spec.title = spec.name;
        }
        specs.push_back(std::move(spec));
    }

    if (specs.empty())
    {
        throw std::runtime_error("read_template_spec_1d_tsv: no templates read from " + path);
    }

    return specs;
}

} // namespace nuxsec
