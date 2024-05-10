#include <fstream>
#include <iomanip>

#include "registry.hpp"

void Registry::gen_c_module(const Module &mod, const std::string &out_dir)
{
    auto file_name = mod.name + "_Types.h";
    auto file_path = out_dir + "/" + file_name;
    std::string indent("\n");

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
    w << "//!STARTOFREGISTRYGENERATEDFILE '" << file_name << "'";
    w << indent << "//!";
    w << indent << "//! WARNING This file is generated automatically by the FAST registry.";
    w << indent << "//! Do not edit.  Your changes to this file will be lost.";
    w << indent << "//!";
    w << indent;
    w << indent << "#ifndef _" << mod.name << "_TYPES_H";
    w << indent << "#define _" << mod.name << "_TYPES_H";
    w << indent;
    w << indent << "#ifdef _WIN32 //define something for Windows (32-bit)";
    w << indent << "\t#include \"stdbool.h\"";
    w << indent << "\t#define CALL __declspec(dllexport)";
    w << indent << "#elif _WIN64 //define something for Windows (64-bit)";
    w << indent << "\t#include \"stdbool.h\"";
    w << indent << "\t#define CALL __declspec(dllexport) ";
    w << indent << "#else";
    w << indent << "\t#include <stdbool.h>";
    w << indent << "\t#define CALL ";
    w << indent << "#endif";

    // Loop through data types in module
    for (auto &dt_name : mod.ddt_names)
    {
        // Get derive data types in module
        auto it = mod.data_types.find(dt_name);
        auto &dt = *it->second;
        if (dt.tag != DataType::Tag::Derived)
            continue;
        auto &ddt = dt.derived;

        // Start of struct
        w << indent;
        w << indent << "typedef struct " << ddt.type_fortran << " {";
        indent += "\t";
        w << indent << "void *object;";

        // Loop through fields
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
                    w << indent << std::setw(28) << std::left << field.data_type->c_type() + " *" + field.name + ";"
                      << "int " << field.name << "_Len;";
                }
                else if (field.data_type->tag == DataType::Tag::Character)
                {
                    if (field.rank == 0)
                    {
                        w << indent << field.data_type->c_type() << " " << field.name << "["
                          << field.data_type->basic.string_len << "];";
                    }
                }
                else
                {
                    w << indent << field.data_type->c_type() << " " << field.name << ";";
                }
            }
            for (int i = 0; i < field.rank; i++)
            {
                if (!field.is_allocatable &&
                    (field.data_type->tag != DataType::Tag::Character || field.rank == 0))
                    w << "[" << field.dims[i].upper_bound << "-" << field.dims[i].lower_bound << "+1];";
            }
        }

        indent.erase(indent.size() - 1);
        w << indent << "} " << ddt.type_fortran << "_t;";
    }

    //--------------------------------------------------------------------------
    // Write struct containing all of the module's derived types
    //--------------------------------------------------------------------------

    w << indent;
    w << indent << "typedef struct " << mod.nickname << "_UserData {";
    indent += "\t";

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
        w << indent << std::setw(28) << std::left << ddt.type_fortran + "_t"
          << " " << mod.nickname << "_" << ddt.interface->name_short << ";";
    }

    indent.erase(indent.size() - 1);
    w << indent << "} " << mod.nickname << "_t;";

    // Write file footer
    w << indent;
    w << indent << "#endif // _" << mod.name << "_TYPES_H";
    w << indent;
    w << indent << "//!ENDOFREGISTRYGENERATEDFILE";
    w << indent;
}