#ifndef REGISTRY_HPP
#define REGISTRY_HPP

#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

std::string tolower(std::string s);

// case-independent (ci) string less_than: returns true if s1 < s2
struct ci_less
{
    // case-independent (ci) compare_less binary function
    struct nocase_compare
    {
        bool operator()(const unsigned char &c1, const unsigned char &c2) const
        {
            return tolower(c1) < tolower(c2);
        }
    };
    bool operator()(const std::string &s1, const std::string &s2) const
    {
        return std::lexicographical_compare(s1.begin(), s1.end(), // source range
                                            s2.begin(), s2.end(), // dest range
                                            nocase_compare());    // comparison
    }
};

enum class Period
{
    None,
    TwoPi,
};

struct Module;
struct DataType;

struct InterfaceData
{
    std::string name;
    std::string name_short;
    bool only_reals;

    InterfaceData(std::string name, std::string name_short, bool only_reals)
        : name(name), name_short(name_short), only_reals(only_reals)
    {
    }
};

struct DimSpec
{
    size_t i = 0;
    bool is_deferred = false;
    bool is_pointer = false;
    std::string lower_bound = "1";
    std::string upper_bound = "-1";

    DimSpec(std::string spec)
    {
        // Get indices of first colon and asterisk
        auto i = spec.find(":");
        auto j = spec.find("*");

        // If colon was found
        if (i != std::string::npos)
        {
            // If colon is the only character, this is a deferred dimension
            this->is_deferred = spec.size() == 1;

            // If colon isn't first, then parse the lower bound, otherwise 1
            this->lower_bound = i > 0 ? spec.substr(0, i) : "1";

            // Parse the upper bound
            this->upper_bound = this->is_deferred ? "-1" : spec.substr(i + 1);
        }
        // If asterisk was found
        else if (j != std::string::npos)
        {
            this->is_deferred = true;
            this->is_pointer = true;
        }
        // Otherwise, spec contains upper bound
        else
        {
            this->lower_bound = "1";
            this->upper_bound = spec;
        }
    }
};

struct Field
{
    std::string name;
    std::shared_ptr<DataType> data_type;
    std::vector<DimSpec> dims;
    std::string init_value = "";
    std::string desc = "-";
    std::string units = "-";
    Period gen_periodic = Period::None;
    int rank = 0;
    bool is_pointer = false;
    bool is_allocatable = false;
    bool is_target = false;

    Field(const std::string &name, std::shared_ptr<DataType> const &type, const std::string &dims,
          const std::string &ctrl, const std::string &init_value, const std::string &desc,
          const std::string &units)
    {
        if (name[0] == '&')
        {
            this->name = name.substr(1);
            this->is_target = true;
            this->is_pointer = true;
            this->is_allocatable = true;
        }
        else if (name[0] == '*')
        {
            this->name = name.substr(1);
            this->is_pointer = true;
            this->is_allocatable = true;
        }
        else
        {
            this->name = name;
        }

        this->data_type = type;

        if (ctrl.compare("2pi") == 0)
        {
            this->gen_periodic = Period::TwoPi;
        }

        if (desc.compare("-") != 0)
        {
            this->desc = desc;
        }

        if (units.compare("-") != 0)
        {
            this->units = units;
        }

        if (dims.compare("-") != 0)
        {
            // Parse dims, throw exception on error
            if (this->parse_dims(dims) != 0)
            {
                throw std::invalid_argument("invalid dimensions: " + dims);
            }

            // Add dimension number
            for (size_t i = 0; i < this->dims.size(); ++i)
            {
                this->dims[i].i = i + 1;
            }

            // Get field rank (number of dimensions)
            this->rank = static_cast<int>(this->dims.size());

            // Field is a pointer if any dim is a pointer
            this->is_pointer |= std::any_of(this->dims.begin(), this->dims.end(),
                                            [](const DimSpec &ds)
                                            { return ds.is_pointer; });

            // Field is allocatable if any dim is deferred
            this->is_allocatable |= std::any_of(this->dims.begin(), this->dims.end(),
                                                [](const DimSpec &ds)
                                                { return ds.is_deferred; });
        }

        // If field is a pointer, initialize to null
        if (this->is_pointer)
        {
            this->init_value = "null()";
        }
        // If field is allocatable, then no initialization
        else if (this->is_allocatable)
        {
            this->init_value = "";
        }
        // If initialization is not empty
        else if (init_value.compare("-") != 0)
        {
            this->init_value = init_value;
            if (tolower(init_value).compare("f") == 0)
            {
                this->init_value = ".false.";
            }
            else if (tolower(init_value).compare("t") == 0)
            {
                this->init_value = ".true.";
            }
        }
    }

    int parse_dims(std::string dim_field)
    {
        // If no dimensions specified
        if (dim_field.size() == 0)
            return 0;

        // Remove leading and trailing braces
        if (dim_field[0] == '{')
            dim_field = dim_field.substr(1);
        if (dim_field.back() == '}')
            dim_field.pop_back();

        // If dim field is only digits, parse number
        if (std::all_of(dim_field.begin(), dim_field.end(), ::isdigit))
        {
            this->dims.push_back(DimSpec(dim_field));
            return 0;
        }

        // If all dims are colons or asterisks, no braces
        if (std::all_of(dim_field.begin(), dim_field.end(), [](char c)
                        { return c == '*'; }) ||
            std::all_of(dim_field.begin(), dim_field.end(), [](char c)
                        { return c == ':'; }))
        {
            for (auto &dim : dim_field)
            {
                this->dims.push_back(DimSpec(std::string(1, dim)));
            }
            return 0;
        }

        // Split by braces
        std::regex split("\\}\\{");
        std::sregex_token_iterator iter(dim_field.begin(), dim_field.end(), split, -1);
        std::sregex_token_iterator re_end;
        for (; iter != re_end; ++iter)
        {
            this->dims.push_back(DimSpec(*iter));
        }

        return 0;
    }
};

struct DataType
{
    enum class Tag
    {
        Integer,
        Real,
        Logical,
        Character,
        Derived,
    };
    Tag tag;

    struct Basic
    {
        std::string name;
        std::string type_fortran;
        std::string string_len;
        int bit_size = 0;
    };
    Basic basic;

    struct Derived
    {
        std::string name;
        std::string name_short;
        std::string type_fortran;
        std::shared_ptr<Module> module;
        std::vector<Field> fields;
        bool contains_mesh = false;
        std::shared_ptr<InterfaceData> interface;
        int max_rank = 0;

        bool only_contains_reals()
        {
            // Loop through fields
            for (const auto &field : this->fields)
            {
                // Switch based on field data type
                switch (field.data_type->tag)
                {

                // Field is a derived type, so check its fields and
                // return false if it doesn't only contain reals
                case Tag::Derived:
                    if (!field.data_type->derived.only_contains_reals())
                        return false;
                    continue;

                // Field is a real, continue
                case Tag::Real:
                    continue;

                // Field is not a real, return false
                case Tag::Character:
                case Tag::Integer:
                case Tag::Logical:
                    return false;
                }
            }

            // Derived data type and all of its fields only contain reals
            return true;
        }
    };
    Derived derived;

    // Constructor for basic type
    DataType(const std::string &name, const std::string &type_fortran, const Tag &type,
             const int bit_size = 0, const std::string &string_len = "")
        : tag(type)
    {
        this->basic.name = name;
        this->basic.type_fortran = type_fortran;
        this->basic.string_len = string_len;
        this->basic.bit_size = bit_size;
    }

    // Constructor for derived type
    DataType(std::shared_ptr<Module> mod, const std::string &name,
             const std::string &name_short = "", const std::string &name_prefixed = "")
        : tag(Tag::Derived)
    {
        this->derived.name = name;
        this->derived.module = mod;
        this->derived.name_short = name_short.empty() ? name : name_short;
        this->derived.type_fortran = name_prefixed.empty() ? name : name_prefixed;
        this->derived.contains_mesh =
            (tolower(name).compare("meshtype") == 0) || (tolower(name).compare("meshmaptype") == 0);
    }

    std::string c_type()
    {
        switch (this->tag)
        {
        case DataType::Tag::Integer:
            return "int";
        case DataType::Tag::Logical:
            return "bool";
        case DataType::Tag::Character:
            return "char";
        case DataType::Tag::Real:
            switch (this->basic.bit_size)
            {
            case 0:
                return "float";
            case 32:
                return "float";
            case 64:
                return "double";
            }
        case DataType::Tag::Derived:
            return "invalid";
        }
        return "invalid";
    }

    std::string c_types_binding()
    {
        switch (this->tag)
        {
        case DataType::Tag::Integer:
            return "INTEGER(KIND=C_INT)";
        case DataType::Tag::Logical:
            return "LOGICAL(KIND=C_BOOL)";
        case DataType::Tag::Character:
            return "CHARACTER(KIND=C_CHAR), DIMENSION(" + this->basic.string_len + ")";
        case DataType::Tag::Real:
            switch (this->basic.bit_size)
            {
            case 0:
                return "REAL(KIND=C_FLOAT)";
            case 32:
                return "REAL(KIND=C_FLOAT)";
            case 64:
                return "REAL(KIND=C_DOUBLE)";
            }
        case DataType::Tag::Derived:
            return "INVALID";
        }
        return "INVALID";
    }
};

struct Parameter
{
    std::string name;
    std::shared_ptr<DataType> type;
    std::string value = "";
    std::string desc = "-";
    std::string units = "-";

    Parameter(const std::string &name, std::shared_ptr<DataType> &type, const std::string &value,
              const std::string &desc, const std::string &units)
    {
        this->name = name;
        this->type = type;
        if (value.compare("-") != 0)
        {
            this->value = value;
        }
        if (desc.compare("-") != 0)
        {
            this->desc = desc;
        }
        if (desc.compare("-") != 0)
        {
            this->desc = desc;
        }
    }
};

struct Module
{
    std::string name;
    std::string nickname;
    std::vector<Parameter> params;
    std::map<std::string, std::shared_ptr<DataType>, ci_less> data_types;
    std::vector<std::string> ddt_names;
    bool is_root = false;

    Module(std::string name, std::string nickname, bool is_root)
        : name(name), nickname(nickname), is_root(is_root)
    {
    }
};

struct Registry
{
    std::vector<std::string> include_dirs = {"."};
    std::set<std::string> include_files;
    std::vector<std::string> use_modules;
    std::map<std::string, std::shared_ptr<InterfaceData>, ci_less> interface_map;
    std::map<std::string, std::shared_ptr<Module>, ci_less> modules;
    std::map<std::string, std::shared_ptr<DataType>, ci_less> data_types;
    bool gen_c_code = false;
    bool no_extrap_interp = false;
    bool gen_inc_subs = false;

    Registry()
    {
        // Basic types
        auto IntKi =
            std::make_shared<DataType>("IntKi", "INTEGER(IntKi)", DataType::Tag::Integer, 32);
        auto SiKi = std::make_shared<DataType>("SiKi", "REAL(SiKi)", DataType::Tag::Real, 32);
        auto R4Ki = std::make_shared<DataType>("R4Ki", "REAL(R4Ki)", DataType::Tag::Real, 32);
        auto ReKi = std::make_shared<DataType>("ReKi", "REAL(ReKi)", DataType::Tag::Real);
        auto R8Ki = std::make_shared<DataType>("R8Ki", "REAL(R8Ki)", DataType::Tag::Real, 64);
        auto DbKi = std::make_shared<DataType>("DbKi", "REAL(DbKi)", DataType::Tag::Real, 64);
        auto logical = std::make_shared<DataType>("Logical", "LOGICAL", DataType::Tag::Logical);

        // Derived types
        auto mesh = std::make_shared<DataType>(nullptr, "MeshType", "MeshType", "MeshType");
        auto dll = std::make_shared<DataType>(nullptr, "DLL_Type");

        // Map of data types
        this->data_types = std::map<std::string, std::shared_ptr<DataType>, ci_less>{
            {"integer", IntKi},
            {"intki", IntKi},
            {"b4ki", IntKi},
            {"real", ReKi},
            {"reki", ReKi},
            {"siki", SiKi},
            {"r4ki", R4Ki},
            {"r8ki", R8Ki},
            {"doubleprecision", DbKi},
            {"dbki", DbKi},
            {"logical", logical},
            {"meshtype", mesh},
            {"dll_type", dll},
        };

        this->interface_map = std::map<std::string, std::shared_ptr<InterfaceData>, ci_less>{
            {"InitInputType", std::make_shared<InterfaceData>("InitInputType", "InitInput", false)},
            {"InitOutputType",
             std::make_shared<InterfaceData>("InitOutputType", "InitOutput", false)},
            {"InputType", std::make_shared<InterfaceData>("InputType", "Input", true)},
            {"OutputType", std::make_shared<InterfaceData>("OutputType", "Output", true)},
            {"ContinuousStateType",
             std::make_shared<InterfaceData>("ContinuousStateType", "ContState", true)},
            {"DiscreteStateType",
             std::make_shared<InterfaceData>("DiscreteStateType", "DiscState", true)},
            {"ConstraintStateType",
             std::make_shared<InterfaceData>("ConstraintStateType", "ConstrState", true)},
            {"OtherStateType",
             std::make_shared<InterfaceData>("OtherStateType", "OtherState", false)},
            {"MiscVarType", std::make_shared<InterfaceData>("MiscVarType", "Misc", false)},
            {"ParameterType", std::make_shared<InterfaceData>("ParameterType", "Param", false)},
            {"PartialOutputPInputType",
             std::make_shared<InterfaceData>("PartialOutputPInputType", "dYdu", true)},
            {"PartialContStatePInputType",
             std::make_shared<InterfaceData>("PartialContStatePInputType", "dXdu", true)},
            {"PartialDiscStatePInputType",
             std::make_shared<InterfaceData>("PartialDiscStatePInputType", "dXddu", true)},
            {"PartialConstrStatePInputType",
             std::make_shared<InterfaceData>("PartialConstrStatePInputType", "dZdu", true)},
        };
    }

    // Parsing
    void parse(const std::string &file_name, const int recurse_level);
    int parse_line(const std::string &line, std::vector<std::string> &fields_prev,
                   const int recurse_level);
    std::shared_ptr<DataType> find_data_type(const std::string &type_name,
                                             std::shared_ptr<Module> mod = nullptr)
    {
        // Pointer to type
        std::shared_ptr<DataType> data_type;

        // Get map of data types to search
        // If module was provided, search it; otherwise, search registry
        auto &data_types = mod == nullptr ? this->data_types : mod->data_types;

        // Search for type in registry, return if found
        auto it = data_types.find(type_name);
        if (it != data_types.end())
        {
            return it->second;
        }

        // If type starts with character (string type), build type and return it
        if (tolower(type_name).compare(0, 9, "character") == 0)
        {
            // Get string length
            auto string_len = type_name.substr(10, type_name.size() - 11);

            // Build type
            data_type = std::make_shared<DataType>(type_name, type_name, DataType::Tag::Character,
                                                   0, string_len);

            // Add type to registry
            this->data_types[type_name] = data_type;
            return data_type;
        }

        return nullptr;
    }

    // Output
    void gen_module_files(std::string const &out_dir);
    void gen_fortran_module(const Module &mod, const std::string &out_dir);
    void gen_c_module(const Module &mod, const std::string &out_dir);
    void gen_fortran_subs(std::ostream &w, const Module &mod);
};

#endif
