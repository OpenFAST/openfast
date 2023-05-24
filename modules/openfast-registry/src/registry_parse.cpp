#include <fstream>
#include <iomanip>

#include "registry.hpp"

const int MAX_FIELDS = 10;

void Registry::parse(const std::string &file_name, const int recurse_level)
{
    std::ifstream inp_file;
    std::vector<std::string> fields_prev;

    // If this is the root file, open given file name
    if (recurse_level == 0)
    {
        std::cerr << "input file: " << file_name << std::endl;
        inp_file.open(file_name);
        if (!inp_file)
        {
            std::cerr << "Registry program cannot open " << file_name << " for reading. ";
            std::cerr << "Ending." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    // Otherwise, find and open include file
    else
    {
        // If this include file has been parsed, return
        if (this->include_files.find(file_name) != this->include_files.end())
            return;

        // Loop through directories and try to open file, break on success
        for (auto &dir : this->include_dirs)
        {
            inp_file.open(dir + "/" + file_name);
            if (inp_file)
                break;
        }

        // If file not opened successfully, exit
        if (!inp_file)
        {
            std::cerr << "Registry error: cannot open '" << file_name << "'." << std::endl;
            exit(EXIT_FAILURE);
        }

        // Display message about opening file
        std::cerr << "opening " << file_name << std::endl;

        // Add file to list of includes
        this->include_files.insert(file_name);
    }

    // Loop through lines in file and parse
    std::string line;
    for (size_t line_num = 1; std::getline(inp_file, line); ++line_num)
    {
        // Parse line into record
        if (this->parse_line(line, fields_prev, recurse_level) != 0)
        {
            std::cerr << "Error reading " << file_name << ":" << line_num << "\n";
            exit(EXIT_FAILURE);
        }
    }

    // If this file is directly included by the root file, save use module
    if (recurse_level == 1)
    {
        auto slash_index = fields_prev[1].find("/");
        bool has_slash = slash_index != std::string::npos;
        auto module_name = has_slash ? fields_prev[1].substr(0, slash_index) : fields_prev[1];
        this->use_modules.push_back(module_name);
    }
}

int Registry::parse_line(const std::string &line, std::vector<std::string> &fields_prev,
                         const int recurse_level)
{
    std::istringstream iss(line);
    std::string s;
    std::vector<std::string> fields;

    // Read fields from line while respecting quotes
    while (iss >> std::quoted(s))
    {
        // If # found in unquoted field, break iteration
        if (s.find("#") != std::string::npos && s.find(" ") == std::string::npos)
            break;

        fields.push_back(s);
    }

    // Skip empty line
    if (fields.size() == 0 || fields[0][0] == '#')
        return EXIT_SUCCESS;

    //--------------------------------------------------------------------------
    // Include Line
    //--------------------------------------------------------------------------

    if (fields.size() == 2 &&
        (tolower(fields[0]).compare("include") == 0 || tolower(fields[0]).compare("usefrom") == 0))
    {
        auto file_name = fields[1];
        this->parse(file_name, recurse_level + 1);
        return EXIT_SUCCESS;
    }

    //--------------------------------------------------------------------------
    // Populate Fields
    //--------------------------------------------------------------------------

    // Resize and fill remaining fields
    fields.resize(MAX_FIELDS, "-");

    // Propagate field values from previous fields if requested
    for (int i = 0; i < MAX_FIELDS; i++)
        if (fields[i].compare("^") == 0)
            fields[i] = fields_prev[i];

    // Update previous fields to current values
    fields_prev = fields;

    //--------------------------------------------------------------------------
    // Get Module
    //--------------------------------------------------------------------------

    // Shared pointer to module
    std::shared_ptr<Module> mod;

    // Is this the root module
    auto is_root = recurse_level == 0;

    // Parse module name and nickname from field
    auto slash_index = fields[1].find("/");
    bool has_slash = slash_index != std::string::npos;
    auto module_name = has_slash ? fields[1].substr(0, slash_index) : fields[1];
    auto module_nickname = has_slash ? fields[1].substr(slash_index + 1) : fields[1];

    // Find module in map or create and add it to map
    auto it = this->modules.find(module_name);
    if (it == this->modules.end())
    {
        mod = std::make_shared<Module>(module_name, module_nickname, is_root);
        this->modules[module_name] = mod;
    }
    else
    {
        mod = it->second;
    }

    //--------------------------------------------------------------------------
    // Parameter Line
    //--------------------------------------------------------------------------

    if (tolower(fields[0]).compare("param") == 0)
    {
        auto name = fields[4];
        auto type = fields[3];
        auto value = fields[6];
        auto desc = fields[8];
        auto units = fields[9];

        // Find parameter type in registry, display message if not found
        auto param_type = this->find_data_type(type);
        if (param_type == nullptr)
        {
            std::cerr << "Registry error: type " << type << " used before defined for " << name
                      << std::endl;
            return EXIT_FAILURE;
        }

        // Add parameter to module
        mod->params.push_back(Parameter(name, param_type, value, desc, units));
        return EXIT_SUCCESS;
    }

    //--------------------------------------------------------------------------
    // Derived Type Line
    //--------------------------------------------------------------------------

    if ((tolower(fields[0]).compare("typedef") == 0) ||
        (tolower(fields[0]).compare("usefrom") == 0))
    {
        auto ddt_name_base = fields[2];
        auto field_type_name = fields[3];
        auto name = fields[4];
        auto dims = fields[5];
        auto init_value = fields[6];
        auto ctrl = fields[7];
        auto desc = fields[8];
        auto units = fields[9];

        // Get derived data type name
        auto ddt_name = ddt_name_base;
        auto ddt_name_short = ddt_name_base;

        // Remove module prefix from name
        std::string prefix = tolower(mod->nickname) + "_";
        if (tolower(ddt_name_short).compare(0, prefix.size(), prefix) == 0)
        {
            ddt_name_short = ddt_name_short.substr(prefix.size());
        }

        // If interface name was found for derived data type, prepend module nickname
        auto it = this->interface_map.find(ddt_name_short);
        auto is_interface_type = it != this->interface_map.end();
        if (is_interface_type)
        {
            ddt_name = mod->nickname + "_" + ddt_name_short;
        }

        // Get data type from module
        auto ddt_dt = this->find_data_type(ddt_name, mod);

        // If struct type not found and module is not root, get from registry
        if (ddt_dt == nullptr && !mod->is_root)
            ddt_dt = this->find_data_type(ddt_name);

        // If derived data type not found, create and add to module or registry
        if (ddt_dt == nullptr)
        {
            // Get short name from interface if this is an interface type
            if (is_interface_type)
                ddt_name_short = it->second->name_short;

            // Create derived data type
            ddt_dt = std::make_shared<DataType>(mod, ddt_name_base, ddt_name_short, ddt_name);

            // Add interface to type if found
            if (is_interface_type)
                ddt_dt->derived.interface = it->second;

            // Add type module if this is root; otherwise, add to registry
            if (is_root)
            {
                mod->data_types[ddt_name] = ddt_dt;
                mod->ddt_names.push_back(ddt_name);
            }
            else
            {
                this->data_types[ddt_name] = ddt_dt;
            }
        }

        // Get field data type from module or registry
        auto field_dt = this->find_data_type(field_type_name, mod);
        if (field_dt == nullptr)
        {
            field_dt = this->find_data_type(field_type_name);
        }
        if (field_dt == nullptr)
        {
            std::cerr << "Error: type " << field_type_name << " used before defined for " << name
                      << std::endl;
            return EXIT_FAILURE;
        }

        // Create field
        Field field(name, field_dt, dims, ctrl, init_value, desc, units);

        // The field is a target pointer if the following is true:
        // - C code will be generated
        // - The field is allocatable
        // - The field is not a derived type
        // - The field name doesn't start with "writeoutput"
        if (this->gen_c_code && field.is_allocatable &&
            (field.data_type->tag != DataType::Tag::Derived) &&
            (tolower(field.name.substr(0, 11)).compare("writeoutput") != 0))
        {
            field.is_pointer = true;
            field.is_target = true;
        }

        // If field is a mesh derived type (MeshType or MeshMapType)
        // or a derived type that contains a mesh,
        // set flag in derived data type
        if ((field.data_type->tag == DataType::Tag::Derived) &&
            field.data_type->derived.contains_mesh)
            ddt_dt->derived.contains_mesh = true;

        // Accumulate max rank of fields in derived data type
        ddt_dt->derived.max_rank = std::max(ddt_dt->derived.max_rank, field.rank);

        // Add field to derived data type
        ddt_dt->derived.fields.push_back(field);
        return EXIT_SUCCESS;
    }

    // Line is invalid
    std::cerr << "Error: invalid line: '" << line << "'\n";
    return EXIT_FAILURE;
}
