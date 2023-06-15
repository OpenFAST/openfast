#include <fstream>
#include <iomanip>

#include "registry.hpp"

void Registry::gen_c_module(const Module &mod, const std::string &out_dir)
{
    auto file_name = mod.name + "_Types.h";
    auto file_path = out_dir + "/" + file_name;

    // Write message that file is being generated
    std::cerr << "generating " << file_name << std::endl;

    // Open output file, return if error
    std::ofstream w(file_path);
    if (!w)
    {
        std::cerr << "Error creating module file: '" << file_path << "'" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Write file header
    w << "//!STARTOFREGISTRYGENERATEDFILE '" << file_name << "'\n";
    w << "//!\n";
    w << "//! WARNING This file is generated automatically by the FAST registry.\n";
    w << "//! Do not edit.  Your changes to this file will be lost.\n";
    w << "//!\n";
    w << "\n";
    w << "#ifndef _" << mod.name << "_TYPES_H\n";
    w << "#define _" << mod.name << "_TYPES_H\n\n";
    w << "\n";
    w << "#ifdef _WIN32 //define something for Windows (32-bit)\n";
    w << "#  include \"stdbool.h\"\n";
    w << "#  define CALL __declspec( dllexport )\n";
    w << "#elif _WIN64 //define something for Windows (64-bit)\n";
    w << "#  include \"stdbool.h\"\n";
    w << "#  define CALL __declspec( dllexport ) \n";
    w << "#else\n";
    w << "#  include <stdbool.h>\n";
    w << "#  define CALL \n";
    w << "#endif\n";
    w << "\n\n";

    // Loop through data types in module
    for (auto &dt_name : mod.ddt_names)
    {
        // Get derive data types in module
        auto it = mod.data_types.find(dt_name);
        auto &dt = *it->second;
        if (dt.tag != DataType::Tag::Derived)
            continue;
        auto &ddt = dt.derived;

        w << "  typedef struct " << ddt.type_fortran << " {\n";
        w << "    void * object ;\n";
        for (const auto &field : ddt.fields)
        {
            if (field.data_type->tag == DataType::Tag::Derived)
            {
                // TODO:Support derived types
            }
            else // Basic Type
            {
                if (field.is_allocatable)
                {
                    w << "    " << field.data_type->c_type() << " * " << field.name << " ;     int "
                      << field.name << "_Len ;";
                }
                else if (field.data_type->tag == DataType::Tag::Character)
                {
                    if (field.rank == 0)
                    {
                        w << "    " << field.data_type->c_type() << " " << field.name << "["
                          << field.data_type->basic.string_len << "] ;";
                    }
                }
                else
                {
                    w << "    " << field.data_type->c_type() << " " << field.name << " ;";
                }
            }
            for (int i = 0; i < field.rank; i++)
            {
                if (!field.is_allocatable &&
                    (field.data_type->tag != DataType::Tag::Character || field.rank == 0))
                    w << "[" << field.dims[i].upper_bound << "-" << field.dims[i].lower_bound
                      << "+1] ;";
            }
            w << "\n";
        }
        w << "  } " << ddt.type_fortran << "_t ;\n";
    }

    // Write struct containing all of the module's derived types
    w << "  typedef struct " << mod.nickname << "_UserData {\n";
    for (auto &dt_name : mod.ddt_names)
    {
        // Get derived data types with interfaces
        auto it = mod.data_types.find(dt_name);
        auto &dt = *it->second;
        if (dt.tag != DataType::Tag::Derived)
            continue;
        auto &ddt = dt.derived;
        if (ddt.interface == nullptr)
            continue;

        // Write name
        w << "    " << std::setw(30) << std::left << ddt.type_fortran + "_t"
          << " " << mod.nickname << "_" << ddt.interface->name_short << " ;\n";
    }
    w << "  } " << mod.nickname << "_t ;\n";

    // Write file footer
    w << "\n#endif // _" << mod.name << "_TYPES_H\n\n\n";
    w << "//!ENDOFREGISTRYGENERATEDFILE\n";
}