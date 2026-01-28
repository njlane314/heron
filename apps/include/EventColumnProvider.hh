/* -- C++ -- */
/**
 *  @file  apps/include/EventColumnProvider.hh
 *
 *  @brief Column provider for event output schemas.
 */

#ifndef NUXSEC_APPS_EVENT_COLUMN_PROVIDER_H
#define NUXSEC_APPS_EVENT_COLUMN_PROVIDER_H

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "AppUtils.hh"

namespace nuxsec
{

namespace app
{

namespace event
{

/** \brief Provide event columns from compiled defaults or a TSV file. */
class EventColumnProvider
{
  public:
    EventColumnProvider(std::vector<std::string> columns,
                        std::vector<std::pair<std::string, std::string>> schema_columns,
                        std::string columns_tsv_path)
        : m_columns(std::move(columns)),
          m_schema_columns(std::move(schema_columns)),
          m_schema_tag(columns_tsv_path.empty() ? "compiled" : "columns")
    {
        if (!columns_tsv_path.empty())
        {
            load_columns_tsv(columns_tsv_path);
        }
    }

    const std::vector<std::string> &columns() const noexcept { return m_columns; }
    const std::vector<std::pair<std::string, std::string>> &schema_columns() const noexcept
    {
        return m_schema_columns;
    }
    const std::string &schema_tag() const noexcept { return m_schema_tag; }

    std::string schema_tsv() const
    {
        if (m_schema_columns.empty())
        {
            return "";
        }

        std::ostringstream schema;
        schema << "type\tname\n";
        for (const auto &entry : m_schema_columns)
        {
            schema << entry.first << "\t" << entry.second << "\n";
        }
        return schema.str();
    }

  private:
    static std::vector<std::string> split_tsv_line(const std::string &line)
    {
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream stream(line);
        while (std::getline(stream, token, '\t'))
        {
            tokens.push_back(token);
        }
        return tokens;
    }

    static std::string lower(std::string value)
    {
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c)
                       { return static_cast<char>(std::tolower(c)); });
        return value;
    }

    void load_columns_tsv(const std::string &path)
    {
        std::ifstream input(path);
        if (!input)
        {
            throw std::runtime_error("EventColumnProvider: failed to open columns TSV: " + path);
        }

        std::vector<std::string> columns;
        std::vector<std::pair<std::string, std::string>> schema_columns;
        bool header_checked = false;

        std::string line;
        while (std::getline(input, line))
        {
            line = nuxsec::app::trim(line);
            if (line.empty() || line[0] == '#')
            {
                continue;
            }

            const std::vector<std::string> tokens = split_tsv_line(line);
            if (!header_checked)
            {
                header_checked = true;
                if (tokens.size() >= 2 &&
                    lower(nuxsec::app::trim(tokens.at(0))) == "type" &&
                    lower(nuxsec::app::trim(tokens.at(1))) == "name")
                {
                    continue;
                }
            }

            std::string type;
            std::string name;
            if (tokens.size() == 1)
            {
                name = nuxsec::app::trim(tokens.at(0));
                type = "auto";
            }
            else
            {
                type = nuxsec::app::trim(tokens.at(0));
                name = nuxsec::app::trim(tokens.at(1));
            }

            if (name.empty())
            {
                throw std::runtime_error("EventColumnProvider: missing column name in TSV: " + path);
            }

            columns.push_back(name);
            schema_columns.emplace_back(type, name);
        }

        if (!columns.empty())
        {
            m_columns = std::move(columns);
            m_schema_columns = std::move(schema_columns);
        }
    }

    std::vector<std::string> m_columns;
    std::vector<std::pair<std::string, std::string>> m_schema_columns;
    std::string m_schema_tag;
};

} // namespace event

} // namespace app

} // namespace nuxsec

#endif // NUXSEC_APPS_EVENT_COLUMN_PROVIDER_H
