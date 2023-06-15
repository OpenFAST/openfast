#include <fstream>
#include <algorithm>

#include "registry.hpp"
#include "templates.hpp"

const int MAXRECURSE = 9;

void gen_ExtrapInterp(std::ostream &w, const Module &mod, std::string type_name_long,
                      std::string type_kind, const bool useModPrefix);
void gen_copy(std::ostream &w, const Module &mod, const DataType::Derived &ddt,
              const bool gen_c_code);
void gen_destroy(std::ostream &out, const Module &mod, const DataType::Derived &ddt,
                 const bool gen_c_code);
void gen_pack(std::ostream &out, const Module &mod, const DataType::Derived &ddt,
              const bool gen_c_code);
void gen_unpack(std::ostream &out, const Module &mod, const DataType::Derived &ddt,
                bool gen_c_code);
void gen_copy_c2f(std::ostream &w, const Module &mod, const DataType::Derived &ddt);
void gen_copy_f2c(std::ostream &w, const Module &mod, const DataType::Derived &ddt);

std::string dimstr(size_t d)
{
    switch (d)
    {
    case 0:
        return "";
    case 1:
        return "(i1)";
    case 2:
        return "(i1,i2)";
    case 3:
        return "(i1,i2,i3)";
    case 4:
        return "(i1,i2,i3,i4)";
    case 5:
        return "(i1,i2,i3,i4,i5)";
    }
    return " REGISTRY ERROR TOO MANY DIMS ";
}

std::string dimstr_c(size_t d)
{
    switch (d)
    {
    case 0:
        return "";
    case 1:
        return "[i1]";
    case 2:
        return "[i2][i1]";
    case 3:
        return "[i3][i2][i1]";
    case 4:
        return "[i4][i3][i2][i1]";
    case 5:
        return "[i5][i4][i3][i2][i1]";
    }
    return " REGISTRY ERROR TOO MANY DIMS ";
}

void Registry::gen_fortran_module(const Module &mod, const std::string &out_dir)
{
    // Create file name and path
    auto file_name = mod.name + "_Types.f90";
    auto file_path = out_dir + "/" + file_name;
    std::cerr << "generating " << file_name << std::endl;
    bool is_NWTC_Library = false;

    // Open file, exit if error
    std::ofstream w(file_path);
    if (!w)
    {
        std::cerr << "Error creating module file: '" << file_path << "'" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Write preamble
    w << std::regex_replace(FAST_preamble, std::regex("ModuleName"), mod.name);

    // Output USE statements for non-root modules
    for (auto const &mod : this->use_modules)
        if (tolower(mod).compare("nwtc_library") != 0)
            w << "USE " << mod << "_Types\n";

    // If this is the NWTC Library, we're not going to print "USE NWTC_Library"
    if (tolower(mod.name).compare("nwtc_library") == 0)
        w << "USE SysSubs\n"
          << "USE ModReg\n";
    else
        w << "USE NWTC_Library\n";

    w << "IMPLICIT NONE\n";

    // Write parameters to file
    for (const auto &param : mod.params)
    {
        w << "    " << param.type->basic.type_fortran << ", PUBLIC, PARAMETER  :: " << param.name;

        if (!param.value.empty())
            w << " = " << param.value;

        if (param.desc.compare("-") != 0 || param.units.compare("-") != 0)
            w << "      ! " << param.desc << " [" << param.units << "]";

        w << "\n";
    }

    // Loop through data type names in module
    for (auto &dt_name : mod.ddt_names)
    {
        // Get derived data type
        auto &ddt = mod.data_types.find(dt_name)->second->derived;

        // If derived data type should only contain reals,
        // verify that it does, otherwise exit with error
        if ((ddt.interface != nullptr) && ddt.interface->only_reals)
            if (!ddt.only_contains_reals())
                exit(EXIT_FAILURE);

        // Write derived type header
        w << "! =========  " << ddt.type_fortran << (this->gen_c_code ? "_C" : "") << "  =======\n";

        // If requested, write C version of derived data type
        if (this->gen_c_code)
        {
            w << "  TYPE, BIND(C) :: " << ddt.type_fortran << "_C\n";
            w << "   TYPE(C_PTR) :: object = C_NULL_PTR\n";

            for (auto &field : ddt.fields)
            {
                if (field.data_type->tag != DataType::Tag::Derived)
                {
                    if (field.rank == 0)
                    {
                        auto c = field.data_type->c_types_binding();
                        w << "    " << field.data_type->c_types_binding()
                          << " :: " << field.name << " \n";
                    }
                    else
                    {
                        if (field.is_allocatable)
                        {
                            w << "    TYPE(C_ptr) :: " << field.name << " = C_NULL_PTR \n";
                            w << "    INTEGER(C_int) :: " << field.name << "_Len = 0 \n";
                        }
                        else if (field.data_type->tag != DataType::Tag::Character)
                        {
                            w << "    TYPE(C_PTR) :: " << field.name << "(";
                            for (int i = 0; i < field.rank; i++)
                            {
                                w << (i > 0 ? "," : "") << field.dims[i].upper_bound;
                            }
                            w << ")\n";
                        }
                    }
                }
            }
            w << "  END TYPE " << ddt.type_fortran << "_C\n";
        }

        // Write Fortran derived data type
        w << "  TYPE, PUBLIC :: " << ddt.type_fortran << "\n";
        if (this->gen_c_code)
            w << "    TYPE( " << ddt.type_fortran << "_C ) :: C_obj\n";

        // Loop through fields
        for (auto &field : ddt.fields)
        {
            if (field.data_type->tag == DataType::Tag::Derived)
            {
                w << "    TYPE(" << field.data_type->derived.type_fortran << ") ";
            }
            else if (this->gen_c_code && field.is_pointer)
            {
                auto c = field.data_type->c_types_binding();
                w << "    " << field.data_type->c_types_binding() << " ";
            }
            else
            {
                w << "    " << field.data_type->basic.type_fortran << " ";
            }

            if (field.rank > 0)
            {
                w << ", DIMENSION(";

                // If field is allocatable
                if (field.is_allocatable)
                {
                    for (int i = 0; i < field.rank; i++)
                        w << (i == 0 ? ":" : ",:");

                    w << "), " << (field.is_pointer ? "POINTER " : "ALLOCATABLE ");
                }
                // Field is not allocatable
                else
                {
                    bool first = true;
                    for (const auto &dim : field.dims)
                    {
                        w << (first ? "" : ",") << dim.lower_bound << ":"
                          << dim.upper_bound;
                        first = false;
                    }
                    w << ") ";
                }
            }
            else if (field.is_pointer)
            {
                w << ", POINTER";
            }

            w << " :: " << field.name << " ";

            if (field.is_pointer)
            {
                w << "=> NULL() ";
            }
            else if (field.rank == 0 && !field.init_value.empty())
            {
                w << "= " << field.init_value << " ";
            }

            if (field.desc.compare("-") != 0 || field.units.compare("-") != 0)
            {
                w << "     !< " << field.desc << " [" << field.units << "]";
            }

            w << "\n";
        }
        w << "  END TYPE " << ddt.type_fortran << "\n";
        w << "! =======================\n";
    }

    w << "CONTAINS\n";

    // Loop through derived data types
    for (auto &dt_name : mod.ddt_names)
    {
        // Get derived data type
        auto &ddt = mod.data_types.find(dt_name)->second->derived;

        // Generate copy, destroy, pack, and unpack routines
        gen_copy(w, mod, ddt, this->gen_c_code);
        gen_destroy(w, mod, ddt, this->gen_c_code);
        gen_pack(w, mod, ddt, this->gen_c_code);
        gen_unpack(w, mod, ddt, this->gen_c_code);

        // If C code generation requested
        if (this->gen_c_code)
        {
            // Generate C <-> Fortran copy functions
            gen_copy_c2f(w, mod, ddt);
            gen_copy_f2c(w, mod, ddt);
        }
    }

    // If this is the AirfoilInfo module, generate routines for Output and UA_BL_Type
    if (tolower(mod.name).compare("airfoilinfo") == 0)
    {
        gen_ExtrapInterp(w, mod, "OutputType", "ReKi", 1);
        gen_ExtrapInterp(w, mod, "UA_BL_Type", "ReKi", 1);
    }
    else if (!this->no_extrap_interp)
    {
        // If this is the DBEMT module make extrap/interp for ElementInput
        if (tolower(mod.name).compare("dbemt") == 0)
            gen_ExtrapInterp(w, mod, "ElementInputType", "DbKi", 1);

        // Generate extrap/interp routines for module input and output types
        gen_ExtrapInterp(w, mod, "InputType", "DbKi", 1);
        gen_ExtrapInterp(w, mod, "OutputType", "DbKi", 1);
    }

    w << "END MODULE " << mod.name << "_Types\n";
    w << "!ENDOFREGISTRYGENERATEDFILE\n";
}

void gen_copy(std::ostream &w, const Module &mod, const DataType::Derived &ddt,
              const bool gen_c_code)
{
    auto routine_name = mod.nickname + "_Copy" + ddt.name_short;
    std::string indent("\n");

    bool has_alloc = std::any_of(ddt.fields.begin(), ddt.fields.end(), [](Field f)
                                 { return f.is_allocatable; });
    bool has_ddt = std::any_of(ddt.fields.begin(), ddt.fields.end(), [](Field f)
                               { return f.data_type->tag == DataType::Tag::Derived; });
    bool has_ddt_arr = std::any_of(ddt.fields.begin(), ddt.fields.end(), [](Field f)
                                   { return f.data_type->tag == DataType::Tag::Derived && f.rank > 0; });

    w << indent << "subroutine " << routine_name << "(Src" << ddt.name_short
      << "Data, Dst" << ddt.name_short << "Data, CtrlCode, ErrStat, ErrMsg)";
    indent += "   ";
    w << indent << "type(" << ddt.type_fortran << "), intent(" << (ddt.contains_mesh ? "inout" : "in")
      << ") :: Src" << ddt.name_short << "Data";
    w << indent << "type(" << ddt.type_fortran << "), intent(inout) :: Dst" << ddt.name_short << "Data";
    w << indent << "integer(IntKi),  intent(in   ) :: CtrlCode";
    w << indent << "integer(IntKi),  intent(  out) :: ErrStat";
    w << indent << "character(*),    intent(  out) :: ErrMsg";
    if (has_ddt_arr)
    {
        w << indent << "integer(IntKi)  :: ";
        for (int i = 1; i <= ddt.max_rank; i++)
            w << (i > 1 ? ", " : "") << "i" << i;
        w << "";
    }
    if (has_ddt_arr || has_alloc)
        w << indent << "integer(IntKi)                 :: LB(" << ddt.max_rank << "), UB(" << ddt.max_rank << ")";
    if (has_ddt || has_alloc)
        w << indent << "integer(IntKi)                 :: ErrStat2";
    if (has_ddt)
        w << indent << "character(ErrMsgLen)           :: ErrMsg2";
    w << indent << "character(*), parameter        :: RoutineName = '" << routine_name << "'";
    w << indent << "ErrStat = ErrID_None";
    w << indent << "ErrMsg  = ''";

    // Loop through fields
    for (auto &field : ddt.fields)
    {
        std::string alloc_assoc = field.is_pointer ? "associated" : "allocated";
        std::string src = "Src" + ddt.name_short + "Data%" + field.name;
        std::string dst = "Dst" + ddt.name_short + "Data%" + field.name;

        // w << indent << "! " << field.name;

        // If field is a non-target pointer, associate the destination
        // pointer with the source pointer
        if (field.is_pointer && !field.is_target)
        {
            w << indent << dst << " => " << src;
            continue;
        }

        // If field is allocatable
        if (field.is_allocatable)
        {
            w << indent << "if (" << alloc_assoc << "(" << src << ")) then";
            indent += "   ";

            std::string dims("");
            if (field.rank > 0)
            {
                w << indent << "LB(1:" << field.rank << ") = lbound(" << src << ")";
                w << indent << "UB(1:" << field.rank << ") = ubound(" << src << ")";
                for (int d = 1; d <= field.rank; d++)
                    dims += ",LB(" + std::to_string(d) + "):UB(" + std::to_string(d) + ")";
                dims = "(" + dims.substr(1) + ")";
            }

            // If dst alloc/assoc
            w << indent << "if (.not. " << alloc_assoc << "(" << dst << ")) then";
            indent += "   ";
            w << indent << "allocate(" << dst << dims << ", stat=ErrStat2)";
            w << indent << "if (ErrStat2 /= 0) then";
            w << indent << "   call SetErrStat(ErrID_Fatal, 'Error allocating " << dst << ".', ErrStat, ErrMsg, RoutineName)";
            w << indent << "   return";
            w << indent << "end if";

            // bjj: this needs to be updated if we've got multidimensional arrays
            if (gen_c_code && field.is_pointer)
            {
                std::string dst_c = "Dst" + ddt.name_short + "Data%C_obj%" + field.name;
                w << indent << dst_c << "_Len = size(" << dst << ")";
                w << indent << "if (" << dst_c << "_Len > 0) &";
                w << indent << "   " << dst_c << " = c_loc(" << dst << "(";
                for (int d = 1; d <= field.rank; d++)
                    w << (d > 1 ? "," : "") << "LB(" << d << ")";
                w << "))";
            }

            // End if dst alloc/assoc
            indent.erase(indent.size() - 3);
            w << indent << "end if";
        }

        // If derived data type (includes mesh and dll_type)
        if (field.data_type->tag == DataType::Tag::Derived)
        {
            auto &ddt = field.data_type->derived;

            // Get bounds for non-allocated field
            if (field.rank > 0 && !field.is_allocatable)
            {
                w << indent << "LB(1:" << field.rank << ") = lbound(" << src << ")";
                w << indent << "UB(1:" << field.rank << ") = ubound(" << src << ")";
            }

            for (int d = field.rank; d >= 1; d--)
            {
                w << indent << "do i" << d << " = LB(" << d << "), UB(" << d << ")";
                indent += "   ";
            }

            if (ddt.name_short.compare("MeshType") == 0)
            {
                w << indent << "call MeshCopy(" << src << dimstr(field.rank) << ", " << dst
                  << dimstr(field.rank) << ", CtrlCode, ErrStat2, ErrMsg2 )";
                w << indent << "call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)";
                w << indent << "if (ErrStat >= AbortErrLev) return";
            }
            else if (ddt.name_short.compare("DLL_Type") == 0)
            {
                w << indent << dst << " = " << src << "";
            }
            else
            {
                w << indent << "call " << ddt.module->nickname << "_Copy" << ddt.name_short << "("
                  << src << dimstr(field.rank) << ", " << dst << dimstr(field.rank)
                  << ", CtrlCode, ErrStat2, ErrMsg2)";
                w << indent << "call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)";
                w << indent << "if (ErrStat >= AbortErrLev) return";
            }

            for (auto &d : field.dims)
            {
                indent.erase(indent.size() - 3);
                w << indent << "end do";
            }
        }
        else
        {
            // Copy values
            w << indent << dst << " = " << src;

            // If C code and field isn't a pointer, copy data to C object
            if (gen_c_code && !field.is_pointer)
            {
                if (field.rank == 0) // scalar of any type OR a character array
                {
                    std::string tmp = ddt.name_short + "Data%C_obj%" + field.name;
                    w << indent << "Dst" << tmp << " = Src" << tmp;
                }
            }
        }

        // End if for source is allocated/associated
        // If source is not allocated/associated, but destination is allocated
        if (field.is_allocatable)
        {
            indent.erase(indent.size() - 3);
            w << indent << "end if";
        }
    }

    indent.erase(indent.size() - 3);
    w << indent << "end subroutine";
    w << indent;
}

void gen_destroy(std::ostream &w, const Module &mod, const DataType::Derived &ddt,
                 const bool gen_c_code)
{
    auto ddt_data = ddt.name_short + "Data";
    auto routine_name = mod.nickname + "_Destroy" + ddt.name_short;
    std::string indent("\n");

    bool has_ddt = std::any_of(ddt.fields.begin(), ddt.fields.end(), [](Field f)
                               { return f.data_type->tag == DataType::Tag::Derived; });
    bool has_ddt_arr = std::any_of(ddt.fields.begin(), ddt.fields.end(), [](Field f)
                                   { return f.data_type->tag == DataType::Tag::Derived && f.rank > 0; });

    w << indent << "subroutine " << routine_name << "(" << ddt_data << ", ErrStat, ErrMsg)";
    indent += "   ";
    w << indent << "type(" << ddt.type_fortran << "), intent(inout) :: " << ddt_data;
    w << indent << "integer(IntKi),  intent(  out) :: ErrStat";
    w << indent << "character(*),    intent(  out) :: ErrMsg";
    if (has_ddt_arr)
    {
        w << indent << "integer(IntKi)  :: ";
        for (int i = 1; i <= ddt.max_rank; i++)
            w << (i > 1 ? ", " : "") << "i" << i;
        w << indent << "integer(IntKi)  :: LB(" << ddt.max_rank << "), UB(" << ddt.max_rank << ")";
    }
    if (has_ddt)
    {
        w << indent << "integer(IntKi)                 :: ErrStat2";
        w << indent << "character(ErrMsgLen)           :: ErrMsg2";
    }
    w << indent << "character(*), parameter        :: RoutineName = '" << routine_name << "'";
    w << indent << "ErrStat = ErrID_None";
    w << indent << "ErrMsg  = ''";

    // Loop through fields in derived data type
    for (auto &field : ddt.fields)
    {
        auto var = ddt_data + "%" + field.name;
        std::string alloc_assoc = field.is_pointer ? "associated" : "allocated";

        // w << indent << "! " << field.name;

        // If non-target pointer field, just nullify pointer
        if (field.is_pointer && !field.is_target)
        {
            w << indent << "nullify(" << var << ")";
            continue;
        }

        // If field is allocatable
        if (field.is_allocatable)
        {
            w << indent << "if (" << alloc_assoc << "(" << var << ")) then";
            indent += "   ";
        }

        // If field is a derived data type, loop through elements and destroy
        if (field.data_type->tag == DataType::Tag::Derived)
        {
            auto var_dims = var + dimstr(field.rank);

            if (field.rank > 0)
            {
                w << indent << "LB(1:" << field.rank << ") = lbound(" << var << ")";
                w << indent << "UB(1:" << field.rank << ") = ubound(" << var << ")";
            }
            for (int d = field.rank; d >= 1; d--)
            {
                w << indent << "do i" << d << " = LB(" << d << "), UB(" << d << ")";
                indent += "   ";
            }

            if (field.data_type->derived.name.compare("MeshType") == 0)
            {
                w << indent << "call MeshDestroy( " << var_dims << ", ErrStat2, ErrMsg2)";
                w << indent << "call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)";
            }
            else if (field.data_type->derived.name.compare("DLL_Type") == 0)
            {
                w << indent << "call FreeDynamicLib( " << var_dims << ", ErrStat2, ErrMsg2)";
                w << indent << "call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)";
            }
            else
            {
                w << indent << "call " << field.data_type->derived.module->nickname << "_Destroy"
                  << field.data_type->derived.name_short << "(" << var_dims << ", ErrStat2, ErrMsg2)";
                w << indent << "call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)";
            }

            // Close for loops
            for (int d = field.rank; d >= 1; d--)
            {
                indent.erase(indent.size() - 3);
                w << indent << "end do";
            }
        }

        if (field.is_allocatable)
        {
            w << indent << "deallocate(" << var << ")";
            if (field.is_pointer)
            {
                w << indent << var << " => null()";
            }

            if (gen_c_code && field.is_pointer)
            {
                auto var_c = ddt_data + "%C_obj%" + field.name;
                w << indent << var_c << " = c_null_ptr";
                w << indent << var_c << "_Len = 0";
            }

            indent.erase(indent.size() - 3);
            w << indent << "end if";
        }
    }

    indent.erase(indent.size() - 3);
    w << indent << "end subroutine";
    w << indent;
}

void gen_pack(std::ostream &w, const Module &mod, const DataType::Derived &ddt,
              const bool gen_c_code)
{
    auto ddt_data = ddt.name_short + "Data";
    auto routine_name = mod.nickname + "_Pack" + ddt.name_short;
    std::string indent("\n");

    bool has_alloc = std::any_of(ddt.fields.begin(), ddt.fields.end(), [](Field f)
                                 { return f.is_allocatable; });
    bool has_ptr = std::any_of(ddt.fields.begin(), ddt.fields.end(), [](Field f)
                               { return f.is_pointer; });
    bool has_ddt_arr = std::any_of(ddt.fields.begin(), ddt.fields.end(), [](Field f)
                                   { return f.data_type->tag == DataType::Tag::Derived && f.rank > 0; });

    w << indent << "subroutine " << routine_name << "(Buf, Indata)";
    indent += "   ";
    w << indent << "type(PackBuffer), intent(inout) :: Buf";
    w << indent << "type(" << ddt.type_fortran << "), intent(in) :: InData";
    w << indent << "character(*), parameter         :: RoutineName = '" << routine_name << "'";
    if (has_ddt_arr)
    {
        w << indent << "integer(IntKi)  :: ";
        for (int i = 1; i <= ddt.max_rank; i++)
            w << (i > 1 ? ", " : "") << "i" << i;
        w << indent << "integer(IntKi)  :: LB(" << ddt.max_rank << "), UB(" << ddt.max_rank << ")";
    }
    if (has_ptr)
    {
        w << indent << "logical         :: PtrInIndex";
    }

    w << indent << "if (Buf%ErrStat >= AbortErrLev) return";

    if (gen_c_code)
    {
        w << indent << "if (c_associated(InData%C_obj%object)) then";
        w << indent << "   call SetErrStat(ErrID_Severe,'C_obj%object cannot be packed.', Buf%ErrStat, Buf%ErrMsg, RoutineName)";
        w << indent << "   return";
        w << indent << "end if";
    }

    // Pack data
    for (auto &field : ddt.fields)
    {
        auto assoc_alloc = field.is_pointer ? "associated" : "allocated";
        auto var = "InData%" + field.name;

        // w << indent << "! " << field.name;

        if (field.is_allocatable)
        {
            w << indent << "call RegPack(Buf, " << assoc_alloc << "(" << var << "))";
            w << indent << "if (" << assoc_alloc << "(" << var << ")) then";
            indent += "   ";
            if (field.rank > 0)
            {
                w << indent << "call RegPackBounds(Buf, " << field.rank << ", lbound(" << var << "), ubound(" << var << "))";
            }
        }
        if (field.is_pointer)
        {
            w << indent << "call RegPackPointer(Buf, c_loc(" << var << "), PtrInIndex)";
            w << indent << "if (.not. PtrInIndex) then";
            indent += "   ";
        }

        //  call individual routines to pack data from subtypes:
        if (field.data_type->tag == DataType::Tag::Derived)
        {
            auto field_dims = var + dimstr(field.rank);

            if (field.rank > 0)
            {
                w << indent << "LB(1:" << field.rank << ") = lbound(" << var << ")";
                w << indent << "UB(1:" << field.rank << ") = ubound(" << var << ")";
            }

            for (int d = field.rank; d >= 1; d--)
            {
                w << indent << "do i" << d << " = LB(" << d << "), UB(" << d << ")";
                indent += "   ";
            }

            if (field.data_type->derived.name.compare("MeshType") == 0)
            {
                w << indent << "call MeshPack(Buf, " << field_dims << ") ";
            }
            else if (field.data_type->derived.name.compare("DLL_Type") == 0)
            {
                w << indent << "call DLLTypePack(Buf, " << field_dims << ") ";
            }
            else
            {
                w << indent << "call " << field.data_type->derived.module->nickname << "_Pack"
                  << field.data_type->derived.name_short << "(Buf, " << field_dims << ") ";
            }

            for (int d = field.rank; d >= 1; d--)
            {
                indent.erase(indent.size() - 3);
                w << indent << "end do";
            }
        }
        else
        {
            // Intrinsic types are handled by generic Pack method on buffer
            w << indent << "call RegPack(Buf, " << var << ")";
        }

        if (field.is_pointer)
        {
            indent.erase(indent.size() - 3);
            w << indent << "end if";
        }

        if (field.is_allocatable)
        {
            indent.erase(indent.size() - 3);
            w << indent << "end if";
        }
    }
    
    // Check for pack errors at end of routine
    w << indent << "if (RegCheckErr(Buf, RoutineName)) return";

    indent.erase(indent.size() - 3);
    w << indent << "end subroutine";
    w << indent;
}

void gen_unpack(std::ostream &w, const Module &mod, const DataType::Derived &ddt, bool gen_c_code)
{
    auto routine_name = mod.nickname + "_UnPack" + ddt.name_short;
    std::string indent("\n");

    bool has_alloc = std::any_of(ddt.fields.begin(), ddt.fields.end(), [](Field f)
                                 { return f.is_allocatable; });
    bool has_ptr = std::any_of(ddt.fields.begin(), ddt.fields.end(), [](Field f)
                               { return f.is_pointer; });
    bool has_ddt_arr = std::any_of(ddt.fields.begin(), ddt.fields.end(), [](Field f)
                                   { return f.data_type->tag == DataType::Tag::Derived && f.rank > 0; });

    w << indent << "subroutine " << routine_name << "(Buf, OutData)";
    indent += "   ";
    w << indent << "type(PackBuffer), intent(inout)    :: Buf";
    w << indent << "type(" << ddt.type_fortran << "), intent(inout) :: OutData";
    w << indent << "character(*), parameter            :: RoutineName = '" << routine_name << "'";
    if (has_ddt_arr)
    {
        w << indent << "integer(IntKi)  :: ";
        for (int i = 1; i <= ddt.max_rank; i++)
            w << (i > 1 ? ", " : "") << "i" << i;
        w << "";
    }
    if (has_ddt_arr || has_alloc)
    {
        w << indent << "integer(IntKi)  :: LB(" << ddt.max_rank << "), UB(" << ddt.max_rank << ")";
    }
    if (has_alloc)
    {
        w << indent << "integer(IntKi)  :: stat";
        w << indent << "logical         :: IsAllocAssoc";
    }
    if (has_ptr)
    {
        w << indent << "integer(IntKi)  :: PtrIdx";
        w << indent << "type(c_ptr)     :: Ptr";
    }
    w << indent << "if (Buf%ErrStat /= ErrID_None) return";

    // BJJ: TODO:  if there are C types, we're going to have to associate with C data structures....

    // Loop through fields and generate code to unpack data
    for (auto &field : ddt.fields)
    {
        auto field_dims = field.name + dimstr(field.rank);
        std::string var = "OutData%" + field.name;
        std::string var_dims = "OutData%" + field.name + dimstr(field.rank);
        std::string var_c = "OutData%C_obj%" + field.name;
        auto assoc_alloc = field.is_pointer ? "associated" : "allocated";

        // w << indent << "! " << field.name << "";

        if (field.is_allocatable)
        {
            w << indent << "if (" << assoc_alloc << "(" << var << ")) deallocate(" << var << ")";
            w << indent << "call RegUnpack(Buf, IsAllocAssoc)";
            w << indent << "if (RegCheckErr(Buf, RoutineName)) return";
            w << indent << "if (IsAllocAssoc) then";
            indent += "   ";
            if (field.rank > 0)
            {
                w << indent << "call RegUnpackBounds(Buf, " << field.rank << ", LB, UB)";
                w << indent << "if (RegCheckErr(Buf, RoutineName)) return";
            }
        }

        if (field.is_pointer)
        {
            w << indent << "call RegUnpackPointer(Buf, Ptr, PtrIdx)";
            w << indent << "if (RegCheckErr(Buf, RoutineName)) return";
            w << indent << "if (c_associated(Ptr)) then";
            if (field.rank == 0)
            {
                w << indent << "   call c_f_pointer(Ptr, " << var << ")";
            }
            else
            {
                auto rank = std::to_string(field.rank);
                w << indent << "   call c_f_pointer(Ptr, " << var << ", UB(1:" << rank << ")-LB(1:" << rank << "))";
                std::string remap_dims;
                for (int d = 1; d <= field.rank; d++)
                    remap_dims += std::string(d > 1 ? "," : "") + "LB(" + std::to_string(d) + "):";
                w << indent << "   " << var << "(" << remap_dims << ") => " << var;
            }
            w << indent << "else";
            indent += "   ";
        }

        if (field.is_allocatable)
        {
            std::string dims;
            for (int d = 1; d <= field.rank; d++)
                dims += std::string(d == 1 ? "(" : "") + "LB(" + std::to_string(d) + ")" +
                        ":UB(" + std::to_string(d) + ")" + (d < field.rank ? "," : ")");
            w << indent << "allocate(" << var << dims << ",stat=stat)";
            w << indent << "if (stat /= 0) then ";
            w << indent << "   call SetErrStat(ErrID_Fatal, 'Error allocating " << var << ".', Buf%ErrStat, Buf%ErrMsg, RoutineName)";
            w << indent << "   return";
            w << indent << "end if";
        }

        // If this is a pointer, set pointer in buffer pointer index
        if (field.is_pointer)
        {
            w << indent << "Buf%Pointers(PtrIdx) = c_loc(" << var << ")";
        }

        // bjj: this needs to be updated if we've got multiple dimension arrays
        if (gen_c_code && field.is_pointer)
        {
            w << indent << var_c << "_Len = size(" << var << ")";
            w << indent << "if (" << var_c << "_Len > 0) " << var_c << " = c_loc(" << var << "(";
            for (int d = 1; d <= field.rank; d++)
                w << (d > 1 ? "," : "") << "LB(" << d << ")";
            w << "))";
        }

        // Call individual routines to unpack data from subtypes:
        if (field.data_type->tag == DataType::Tag::Derived)
        {
            // Get bounds for non-allocated field
            if (field.rank > 0 && !field.is_allocatable)
            {
                w << indent << "LB(1:" << field.rank << ") = lbound(" << var << ")";
                w << indent << "UB(1:" << field.rank << ") = ubound(" << var << ")";
            }

            for (int d = field.rank; d >= 1; d--)
            {
                w << indent << "do i" << d << " = LB(" << d << "), UB(" << d << ")";
                indent += "   ";
            }

            if (field.data_type->derived.name.compare("MeshType") == 0)
            {
                w << indent << "call MeshUnpack(Buf, " << var_dims << ") ! " << field.name << " ";
            }
            else if (field.data_type->derived.name.compare("DLL_Type") == 0)
            {
                w << indent << "call DLLTypeUnpack(Buf, " << var_dims << ") ! " << field.name << " ";
            }
            else
            {
                w << indent << "call " << field.data_type->derived.module->nickname << "_Unpack"
                  << field.data_type->derived.name_short << "(Buf, " << var_dims << ") ! " << field.name << " ";
            }

            for (int d = field.rank; d >= 1; d--)
            {
                indent.erase(indent.size() - 3);
                w << indent << "end do";
            }
        }
        else
        {
            // Intrinsic types are handled by generic unpack method on buffer
            w << indent << "call RegUnpack(Buf, " << var << ")";
            w << indent << "if (RegCheckErr(Buf, RoutineName)) return";

            // need to move scalars and strings to the %c_obj% type, too!
            // compare with copy routine
            if (gen_c_code && !field.is_pointer && field.rank == 0)
            {
                switch (field.data_type->tag)
                {
                case DataType::Tag::Real:
                case DataType::Tag::Integer:
                case DataType::Tag::Logical:
                    w << indent << var_c << " = " << var << "";
                    break;
                case DataType::Tag::Character:
                    w << indent << var_c << " = transfer(" << var << ", " << var_c << " )";
                    break;
                case DataType::Tag::Derived:
                    break;
                }
            }
        }

        if (field.is_pointer)
        {
            indent.erase(indent.size() - 3);
            w << indent << "end if";
        }

        if (field.is_allocatable)
        {
            indent.erase(indent.size() - 3);
            if (field.is_pointer)
            {
                w << indent << "else";
                w << indent << "   " << var << " => null()";
            }
            w << indent << "end if";
        }
    }

    indent.erase(indent.size() - 3);
    w << indent << "end subroutine";
    w << indent;
}

void gen_extint_order(std::ostream &w, const Module &mod, std::string uy, const int order,
                      const Field &field, const std::string &deref, int recurse_level)
{
    std::string indent, tmp;

    if (recurse_level > MAXRECURSE)
    {
        std::cerr << "REGISTRY ERROR: too many levels of array subtypes\n";
        exit(EXIT_FAILURE);
    }

    auto assoc_alloc = field.is_pointer ? "ASSOCIATED" : "ALLOCATED";

    std::string dims = dimstr(field.rank);
    std::string v1 = uy + "1" + deref + "%" + field.name;
    std::string v2 = uy + "2" + deref + "%" + field.name;
    std::string v3 = uy + "3" + deref + "%" + field.name;
    std::string vout = uy + "_out" + deref + "%" + field.name;

    // check if this is an allocatable array:
    if (field.is_allocatable)
    {
        w << "IF (" << assoc_alloc << "(" << vout << ") .AND. " << assoc_alloc << "(" << v1
          << ")) THEN\n";
    }

    if (field.data_type->tag == DataType::Tag::Derived)
    {
        auto &ddt = field.data_type->derived;

        // If this is a type within this module
        if ((ddt.module != nullptr) && (ddt.module->name == mod.name))
        {
            for (auto &sub_field : ddt.fields)
            {
                std::string field_var = deref + "%" + field.name;

                for (int j = field.rank; j > 0; j--)
                {
                    w << "  DO i" << recurse_level << j << " = LBOUND(" << uy << "_out" << field_var
                      << "," << j << "),UBOUND(" << uy << "_out" << field_var << "," << j << ")\n";
                }

                if (field.rank > 0)
                {
                    field_var += "(";
                    for (int j = 1; j <= field.rank; j++)
                    {
                        field_var += "i" + std::to_string(recurse_level) + std::to_string(j);
                        if (j < field.rank)
                        {
                            field_var += ",";
                        }
                    }
                    field_var += ")";
                }

                gen_extint_order(w, mod, uy, order, sub_field, field_var, recurse_level + 1);
                for (int j = field.rank; j > 0; j--)
                {
                    w << "  ENDDO\n";
                }
            }
        }
        else
        {
            for (int j = field.rank; j > 0; j--)
            {
                w << "  DO i" << j << " = LBOUND(" << vout << "," << j << "),UBOUND(" << vout << ","
                  << j << ")\n";
            }

            if (field.data_type->derived.name.compare("MeshType") == 0)
            {
                if (order == 0)
                {
                    w << "      CALL MeshCopy(" << v1 + dims << ", " << vout + dims
                      << ", MESH_UPDATECOPY, ErrStat2, ErrMsg2 )\n";
                }
                else if (order == 1)
                {
                    w << "      CALL MeshExtrapInterp1(" << v1 + dims << ", " << v2 + dims
                      << ", tin, " << vout + dims << ", tin_out, ErrStat2, ErrMsg2 )\n";
                }
                else if (order == 2)
                {
                    w << "      CALL MeshExtrapInterp2(" << v1 + dims << ", " << v2 + dims << ", "
                      << v3 + dims << ", tin, " << vout + dims
                      << ", tin_out, ErrStat2, ErrMsg2 )\n";
                }
            }
            else
            {
                if (order == 0)
                {
                    w << "      CALL " << field.data_type->derived.module->nickname << "_Copy"
                      << field.data_type->derived.name_short << "(" << v1 + dims << ", "
                      << vout + dims << ", MESH_UPDATECOPY, ErrStat2, ErrMsg2 )\n";
                }
                else if (order == 1)
                {
                    w << "      CALL " << field.data_type->derived.module->nickname << "_"
                      << field.data_type->derived.name_short << "_ExtrapInterp1( " << v1 + dims
                      << ", " << v2 + dims << ", tin, " << vout + dims
                      << ", tin_out, ErrStat2, ErrMsg2 )\n";
                }
                else if (order == 2)
                {
                    w << "      CALL " << field.data_type->derived.module->nickname << "_"
                      << field.data_type->derived.name_short << "_ExtrapInterp2( " << v1 + dims
                      << ", " << v2 + dims << ", " << v3 + dims << ", tin, " << vout + dims
                      << ", tin_out, ErrStat2, ErrMsg2 )\n";
                }
            }

            w << "        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)\n";

            for (int j = field.rank; j >= 1; j--)
            {
                w << "   ENDDO\n";
            }
        }
    }
    else if (field.data_type->tag == DataType::Tag::Real)
    {
        if (order == 0)
        {
            // bjj: this should probably have some "IF ALLOCATED" statements around it, but we're
            // just calling the copy routine
            w << "  " << vout << " = " << v1 << "\n";
        }
        else
        {
            indent = "";
        }

        for (int j = field.rank; j > 0; j--)
        {
            w << indent << "  DO i" << j << " = LBOUND(" << vout << "," << j << "),UBOUND(" << vout
              << "," << j << ")\n";
            indent += "  ";
        }

        if (order == 1)
        {
            if (field.gen_periodic == Period::TwoPi)
            {
                w << indent << "  CALL Angles_ExtrapInterp( " << v1 + dims << ", " << v2 + dims
                  << ", tin, " << vout + dims << ", tin_out )\n";
            }
            else
            {
                w << indent << "  b = -(" << v1 + dims << " - " << v2 + dims << ")\n";
                w << indent << "  " << vout + dims << " = " << v1 + dims << " + b * ScaleFactor\n";
            };
        }
        if (order == 2)
        {
            if (field.gen_periodic == Period::TwoPi)
            {
                w << indent << "  CALL Angles_ExtrapInterp( " << v1 + dims << ", " << v2 + dims
                  << ", " << v3 + dims << ", tin, " << vout + dims << ", tin_out )\n";
            }
            else
            {
                w << indent << "  b = (t(3)**2*(" << v1 + dims << " - " << v2 + dims
                  << ") + t(2)**2*(-" << v1 + dims << " + " << v3 + dims << "))* scaleFactor\n ";
                w << indent << " c = ( (t(2)-t(3))*" << v1 + dims << " + t(3)*" << v2 + dims
                  << " - t(2)*" << v3 + dims << " ) * scaleFactor\n";
                w << indent << "  " << vout + dims << " = " << v1 + dims << " + b  + c * t_out\n";
            }
        }
        for (int j = field.rank; j >= 1; j--)
        {
            indent = "";
            for (int i = 1; i < j; i++)
                indent += "  ";
            w << indent << "  END DO\n";
        }
    }

    // check if this is an allocatable array:
    if (field.is_allocatable)
        w << "END IF ! check if allocated\n";
}

void calc_extint_order(std::ostream &w, const Module &mod, const Field &field, int recurse_level,
                       int &max_rank, int &max_nrecurs, int &max_alloc_ndims)
{
    // bjj: make sure this is consistent with logic of gen_extint_order

    // If recursion level is greater than limit, exit with error
    if (recurse_level > MAXRECURSE)
    {
        std::cerr << "REGISTRY ERROR: too many levels of array subtypes\n";
        exit(EXIT_FAILURE);
    }

    // Update max dims based on field rank
    max_rank = std::max(max_rank, field.rank);

    // Switch based on field data type
    switch (field.data_type->tag)
    {
    case DataType::Tag::Derived:

        // If this derived type belongs to this module
        if (field.data_type->derived.module != nullptr &&
            field.data_type->derived.module->name.compare(mod.name) == 0)
        {
            // Update recursion level
            max_nrecurs = std::max(max_nrecurs, recurse_level);

            // Loop through fields and calculate order
            for (const auto &sub_field : field.data_type->derived.fields)
                calc_extint_order(w, mod, sub_field, recurse_level + 1, max_rank, max_nrecurs,
                                  max_alloc_ndims);
        }
        break;

    case DataType::Tag::Real:
        max_alloc_ndims = std::max(max_alloc_ndims, field.rank);
        break;

    default:
        // TODO: handle other field types
        break;
    }
}

void gen_ExtrapInterp1(std::ostream &w, const Module &mod, const DataType::Derived &ddt,
                       std::string &type_kind, std::string &uy, std::string &mod_prefix,
                       const int max_rank, const int max_nrecurs, const int max_alloc_ndims)
{
    w << "\n";
    w << " SUBROUTINE " << mod.nickname << "_" << ddt.name_short << "_ExtrapInterp1(" << uy << "1, "
      << uy << "2, tin, " << uy << "_out, tin_out, ErrStat, ErrMsg )\n";
    w << "!\n";
    w << "! This subroutine calculates a extrapolated (or interpolated) " << ddt.name_short << " "
      << uy << "_out at time t_out, from previous/future time\n";
    w << "! values of " << uy
      << " (which has values associated with times in t).  Order of the interpolation is 1.\n";
    w << "!\n";
    w << "!  f(t) = a + b * t, or\n";
    w << "!\n";
    w << "!  where a and b are determined as the solution to\n";
    w << "!  f(t1) = " << uy << "1, f(t2) = " << uy << "2\n";
    w << "!\n";
    w << "!" << std::string(130, '.') << "\n";
    w << "\n";

    w << " TYPE(" << ddt.type_fortran << "), INTENT(" << (ddt.contains_mesh == 1 ? "INOUT" : "IN")
      << ")  :: " << uy << "1    ! " << ddt.name_short << " at t1 > t2\n";
    w << " TYPE(" << ddt.type_fortran << "), INTENT(" << (ddt.contains_mesh == 1 ? "INOUT" : "IN")
      << ")  :: " << uy << "2    ! " << ddt.name_short << " at t2 \n";
    w << " REAL(" << type_kind
      << "),         INTENT(IN   )          :: tin(2)   ! Times associated with the "
      << ddt.name_short << "s\n";
    w << " TYPE(" << ddt.type_fortran << "), INTENT(INOUT)  :: " << uy << "_out ! "
      << ddt.name_short << " at tin_out\n";
    w << " REAL(" << type_kind
      << "),         INTENT(IN   )          :: tin_out  ! time to be extrap/interp'd to\n";
    w << " INTEGER(IntKi),     INTENT(  OUT)          :: ErrStat  ! Error status of the operation\n";
    w << " CHARACTER(*),       INTENT(  OUT)          :: ErrMsg   ! Error message if ErrStat /= ErrID_None\n";
    w << "   ! local variables\n";
    w << " REAL(" << type_kind
      << ")                                 :: t(2)     ! Times associated with the "
      << ddt.name_short << "s\n";
    w << " REAL(" << type_kind
      << ")                                 :: t_out    ! Time to which to be extrap/interpd\n";
    w << " CHARACTER(*),                    PARAMETER :: RoutineName = '" << mod.nickname << "_"
      << ddt.name_short << "_ExtrapInterp1'\n";

    w << " REAL(DbKi)                                 :: b        ! temporary for extrapolation/interpolation\n";
    w << " REAL(DbKi)                                 :: ScaleFactor ! temporary for extrapolation/interpolation\n";
    w << " INTEGER(IntKi)                             :: ErrStat2 ! local errors\n";
    w << " CHARACTER(ErrMsgLen)                       :: ErrMsg2  ! local errors\n";
    for (int j = 1; j <= max_rank; j++)
    {
        for (int i = 0; i <= max_nrecurs; i++)
        {
            w << " INTEGER                                    :: i" << i << j << "    ! dim" << j
              << " level " << i << " counter variable for arrays of ddts\n";
        }
    }
    for (int j = 1; j <= max_rank; j++)
    {
        w << " INTEGER                                    :: i" << j << "    ! dim" << j
          << " counter variable for arrays\n";
    }
    w << "    ! Initialize ErrStat\n";
    w << " ErrStat = ErrID_None\n";
    w << " ErrMsg  = \"\"\n";
    w << "    ! we'll subtract a constant from the times to resolve some \n";
    w << "    ! numerical issues when t gets large (and to simplify the equations)\n";
    w << " t = tin - tin(1)\n";
    w << " t_out = tin_out - tin(1)\n";
    w << "\n";

    w << "   IF ( EqualRealNos( t(1), t(2) ) ) THEN\n";
    w << "     CALL SetErrStat(ErrID_Fatal, 't(1) must not equal t(2) to avoid a division-by-zero error.', ErrStat, ErrMsg,RoutineName)\n";
    w << "     RETURN\n";
    w << "   END IF\n\n";

    w << "   ScaleFactor = t_out / t(2)" << std::endl;

    for (const auto &field : ddt.fields)
        gen_extint_order(w, mod, uy, 1, field, "", 0);

    w << " END SUBROUTINE " << mod.nickname << "_" << ddt.name_short << "_ExtrapInterp1\n";
    w << "\n";
}

void gen_ExtrapInterp2(std::ostream &w, const Module &mod, const DataType::Derived &ddt,
                       std::string &type_kind, std::string &uy, std::string &modPrefix,
                       const int max_rank, const int max_nrecurs, const int max_alloc_ndims)
{
    w << "\n";
    w << " SUBROUTINE " << mod.nickname << "_" << ddt.name_short << "_ExtrapInterp2(" << uy << "1, "
      << uy << "2, " << uy << "3, tin, " << uy << "_out, tin_out, ErrStat, ErrMsg )\n";
    w << "!\n";
    w << "! This subroutine calculates a extrapolated (or interpolated) " << ddt.name_short << " "
      << uy << "_out at time t_out, from previous/future time\n";
    w << "! values of " << uy
      << " (which has values associated with times in t).  Order of the interpolation is 2.\n";
    w << "!\n";
    w << "!  expressions below based on either\n";
    w << "!\n";
    w << "!  f(t) = a + b * t + c * t**2\n";
    w << "!\n";
    w << "!  where a, b and c are determined as the solution to\n";
    w << "!  f(t1) = " << uy << "1, f(t2) = " << uy << "2, f(t3) = " << uy << "3\n";
    w << "!\n";
    w << "!" << std::string(130, '.') << "\n";
    w << "\n";

    w << " TYPE(" << ddt.type_fortran << "), INTENT(" << (ddt.contains_mesh == 1 ? "INOUT" : "IN")
      << ")  :: " << uy << "1      ! " << ddt.name_short << " at t1 > t2 > t3\n";
    w << " TYPE(" << ddt.type_fortran << "), INTENT(" << (ddt.contains_mesh == 1 ? "INOUT" : "IN")
      << ")  :: " << uy << "2      ! " << ddt.name_short << " at t2 > t3\n";
    w << " TYPE(" << ddt.type_fortran << "), INTENT(" << (ddt.contains_mesh == 1 ? "INOUT" : "IN")
      << ")  :: " << uy << "3      ! " << ddt.name_short << " at t3\n";
    w << " REAL(" << type_kind
      << "),                 INTENT(IN   )  :: tin(3)    ! Times associated with the "
      << ddt.name_short << "s\n";
    w << " TYPE(" << ddt.type_fortran << "), INTENT(INOUT)  :: " << uy << "_out     ! "
      << ddt.name_short << " at tin_out\n";
    w << " REAL(" << type_kind
      << "),                 INTENT(IN   )  :: tin_out   ! time to be extrap/interp'd to\n";

    w << " INTEGER(IntKi),             INTENT(  OUT)  :: ErrStat   ! Error status of the operation\n";
    w << " CHARACTER(*),               INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None\n";
    w << "   ! local variables\n";
    w << " REAL(" << type_kind
      << ")                                 :: t(3)      ! Times associated with the "
      << ddt.name_short << "s\n";
    w << " REAL(" << type_kind
      << ")                                 :: t_out     ! Time to which to be extrap/interpd\n";
    w << " INTEGER(IntKi)                             :: order     ! order of polynomial fit (max 2)\n";

    w << " REAL(DbKi)                                 :: b        ! temporary for extrapolation/interpolation\n";
    w << " REAL(DbKi)                                 :: c        ! temporary for extrapolation/interpolation\n";
    w << " REAL(DbKi)                                 :: ScaleFactor ! temporary for extrapolation/interpolation\n";
    w << " INTEGER(IntKi)                             :: ErrStat2 ! local errors\n";
    w << " CHARACTER(ErrMsgLen)                       :: ErrMsg2  ! local errors\n";
    w << " CHARACTER(*),            PARAMETER         :: RoutineName = '" << mod.nickname << "_"
      << ddt.name_short << "_ExtrapInterp2'\n";
    for (int j = 1; j <= max_rank; j++)
    {
        for (int i = 0; i <= max_nrecurs; i++)
        {
            w << " INTEGER                                    :: i" << i << j << "    ! dim" << j
              << " level " << i << " counter variable for arrays of ddts\n";
        }
    }
    for (int j = 1; j <= max_rank; j++)
    {
        w << " INTEGER                                    :: i" << j << "    ! dim" << j
          << " counter variable for arrays\n";
    }
    w << "    ! Initialize ErrStat\n";
    w << " ErrStat = ErrID_None\n";
    w << " ErrMsg  = \"\"\n";
    w << "    ! we'll subtract a constant from the times to resolve some \n";
    w << "    ! numerical issues when t gets large (and to simplify the equations)\n";
    w << " t = tin - tin(1)\n";
    w << " t_out = tin_out - tin(1)\n";
    w << "\n";

    w << "   IF ( EqualRealNos( t(1), t(2) ) ) THEN\n";
    w << "     CALL SetErrStat(ErrID_Fatal, 't(1) must not equal t(2) to avoid a division-by-zero error.', ErrStat, ErrMsg,RoutineName)\n";
    w << "     RETURN\n";
    w << "   ELSE IF ( EqualRealNos( t(2), t(3) ) ) THEN\n";
    w << "     CALL SetErrStat(ErrID_Fatal, 't(2) must not equal t(3) to avoid a division-by-zero error.', ErrStat, ErrMsg,RoutineName)\n";
    w << "     RETURN\n";
    w << "   ELSE IF ( EqualRealNos( t(1), t(3) ) ) THEN\n";
    w << "     CALL SetErrStat(ErrID_Fatal, 't(1) must not equal t(3) to avoid a division-by-zero error.', ErrStat, ErrMsg,RoutineName)\n";
    w << "     RETURN\n";
    w << "   END IF\n\n";

    w << "   ScaleFactor = t_out / (t(2) * t(3) * (t(2) - t(3)))\n";

    // recursive
    for (const auto &field : ddt.fields)
    {
        gen_extint_order(w, mod, uy, 2, field, "", 0);
    }

    w << " END SUBROUTINE " << mod.nickname << "_" << ddt.name_short << "_ExtrapInterp2\n";
    w << "\n";
}

void gen_ExtrapInterp(std::ostream &w, const Module &mod, std::string type_name_long,
                      std::string type_kind, const bool useModPrefix)
{
    // Get derived data type from module
    std::string modPrefix = useModPrefix ? mod.nickname + "_" : "";
    auto iter = mod.data_types.find(modPrefix + type_name_long);
    if (iter == mod.data_types.end())
        return;
    const auto &dt = iter->second;
    if (dt == nullptr)
        return;
    const auto &ddt = dt->derived;

    std::string uy = tolower(ddt.name_short).compare("output") == 0 ? "y" : "u";

    w << "\n";
    w << " SUBROUTINE " << mod.nickname << "_" << ddt.name_short << "_ExtrapInterp(" << uy
      << ", t, " << uy << "_out, t_out, ErrStat, ErrMsg )\n";
    w << "!\n";
    w << "! This subroutine calculates a extrapolated (or interpolated) " << ddt.name_short << " "
      << uy << "_out at time t_out, from previous/future time\n";
    w << "! values of " << uy
      << " (which has values associated with times in t).  Order of the interpolation is given by the size of "
      << uy << "\n";
    w << "!\n";
    w << "!  expressions below based on either\n";
    w << "!\n";
    w << "!  f(t) = a\n";
    w << "!  f(t) = a + b * t, or\n";
    w << "!  f(t) = a + b * t + c * t**2\n";
    w << "!\n";
    w << "!  where a, b and c are determined as the solution to\n";
    w << "!  f(t1) = " << uy << "1, f(t2) = " << uy << "2, f(t3) = " << uy
      << "3  (as appropriate)\n";
    w << "!\n";
    w << "!" << std::string(130, '.') << "\n";
    w << "\n";
    w << " TYPE(" << ddt.type_fortran << "), INTENT(" << (ddt.contains_mesh == 1 ? "INOUT" : "IN")
      << ")  :: " << uy << "(:) ! " << ddt.name_short << " at t1 > t2 > t3\n";
    w << " REAL(" << type_kind
      << "),                 INTENT(IN   )  :: t(:)           ! Times associated with the "
      << ddt.name_short << "s\n";
    // Intent must be (INOUT) to prevent ALLOCATABLE array arguments in the DDT from
    // being deallocated in this call. See Sec. 5.1.2.7 of Fortran 2003 standard
    w << " TYPE(" << ddt.type_fortran << "), INTENT(INOUT)  :: " << uy << "_out ! "
      << ddt.name_short << " at tin_out\n";
    w << " REAL(" << type_kind
      << "),                 INTENT(IN   )  :: t_out           ! time to be extrap/interp'd to\n";
    w << " INTEGER(IntKi),             INTENT(  OUT)  :: ErrStat         ! Error status of the operation\n";
    w << " CHARACTER(*),               INTENT(  OUT)  :: ErrMsg          ! Error message if ErrStat /= ErrID_None\n";
    w << "   ! local variables\n";
    w << " INTEGER(IntKi)                             :: order           ! order of polynomial fit (max 2)\n";
    w << " INTEGER(IntKi)                             :: ErrStat2        ! local errors\n";
    w << " CHARACTER(ErrMsgLen)                       :: ErrMsg2         ! local errors\n";
    w << " CHARACTER(*),    PARAMETER                 :: RoutineName = '" << mod.nickname << "_"
      << ddt.name_short << "_ExtrapInterp'\n";
    w << "    ! Initialize ErrStat\n";
    w << " ErrStat = ErrID_None\n";
    w << " ErrMsg  = \"\"\n";
    w << " if ( size(t) .ne. size(" << uy << ")) then\n";
    w << "    CALL SetErrStat(ErrID_Fatal,'size(t) must equal size(" << uy
      << ")',ErrStat,ErrMsg,RoutineName)\n";
    w << "    RETURN\n";
    w << " endif\n";
    w << " order = SIZE(" << uy << ") - 1\n";
    w << " IF ( order .eq. 0 ) THEN\n";
    w << "   CALL " << mod.nickname << "_Copy" << ddt.name_short << "(" << uy << "(1), " << uy
      << "_out, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )\n";
    w << "     CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)\n";
    w << " ELSE IF ( order .eq. 1 ) THEN\n";
    w << "   CALL " << mod.nickname << "_" << ddt.name_short << "_ExtrapInterp1(" << uy << "(1), "
      << uy << "(2), t, " << uy << "_out, t_out, ErrStat2, ErrMsg2 )\n";
    w << "     CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)\n";
    w << " ELSE IF ( order .eq. 2 ) THEN\n";
    w << "   CALL " << mod.nickname << "_" << ddt.name_short << "_ExtrapInterp2(" << uy << "(1), "
      << uy << "(2), " << uy << "(3), t, " << uy << "_out, t_out, ErrStat2, ErrMsg2 )\n";
    w << "     CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)\n";
    w << " ELSE \n";
    w << "   CALL SetErrStat(ErrID_Fatal,'size(" << uy
      << ") must be less than 4 (order must be less than 3).',ErrStat,ErrMsg,RoutineName)\n";
    w << "   RETURN\n";
    w << " ENDIF \n";
    w << " END SUBROUTINE " << mod.nickname << "_" << ddt.name_short << "_ExtrapInterp\n";
    w << "\n";

    // bjj: this is max for module, not for type_name_long
    int max_rank = 0;    // mod.module_ddt_list->max_ndims;
    int max_nrecurs = 0; // MAXRECURSE;
    int max_alloc_ndims = 0;

    // Recursively calculate extrap/interp order
    for (const auto &field : ddt.fields)
    {
        calc_extint_order(w, mod, field, 0, max_rank, max_nrecurs, max_alloc_ndims);
    }

    // Generate first order extrap/interp routine
    gen_ExtrapInterp1(w, mod, ddt, type_kind, uy, modPrefix, max_rank, max_nrecurs,
                      max_alloc_ndims);

    // Generate second order extrap/interp routine
    gen_ExtrapInterp2(w, mod, ddt, type_kind, uy, modPrefix, max_rank, max_nrecurs,
                      max_alloc_ndims);
}

void gen_copy_c2f(std::ostream &w, const Module &mod, const DataType::Derived &ddt)
{
    std::string routine_name = mod.nickname + "_C2Fary_Copy" + ddt.name_short;

    w << " SUBROUTINE " << routine_name << "( " << ddt.name_short
      << "Data, ErrStat, ErrMsg, SkipPointers )\n";
    w << "    TYPE(" << ddt.type_fortran << "), INTENT(INOUT) :: " << ddt.name_short << "Data\n";
    w << "    INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n";
    w << "    CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n";
    w << "    LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers\n";
    w << "    ! \n";
    w << "    LOGICAL                        :: SkipPointers_local\n";
    w << "    ErrStat = ErrID_None\n";
    w << "    ErrMsg  = \"\"\n\n";
    w << "    IF (PRESENT(SkipPointers)) THEN\n";
    w << "       SkipPointers_local = SkipPointers\n";
    w << "    ELSE\n";
    w << "       SkipPointers_local = .false.\n";
    w << "    END IF\n";

    // Loop through fields in derived data type
    for (const auto &field : ddt.fields)
    {
        // If field doesn't have a data type, continue
        if (field.data_type == nullptr)
        {
            continue;
        }

        // If field is a derived type, print warning and continue
        if (field.data_type->tag == DataType::Tag::Derived)
        {
            std::cerr << "Registry WARNING: derived data type " << field.name << " of type "
                      << field.data_type->derived.name << " is not passed through C interface\n";
            continue;
        }

        std::string var_f = ddt.name_short + "Data%" + field.name;
        std::string var_c = ddt.name_short + "Data%C_obj%" + field.name;
        if (field.is_pointer)
        {
            w << "\n    ! -- " << field.name << " " << ddt.name_short << " Data fields\n";
            w << "    IF ( .NOT. SkipPointers_local ) THEN\n";
            w << "       IF ( .NOT. C_ASSOCIATED( " << var_c << " ) ) THEN\n";
            w << "          NULLIFY( " << var_f << " )\n";
            w << "       ELSE\n";
            w << "          CALL C_F_POINTER(" << var_c << ", " << var_f << ", (/" << var_c
              << "_Len/))\n";
            w << "       END IF\n";
            w << "    END IF\n";
        }
        else if (!field.is_allocatable)
        {
            switch (field.data_type->tag)
            {
            case DataType::Tag::Real:
            case DataType::Tag::Integer:
            case DataType::Tag::Logical:
                w << "    " << var_f << " = " << var_c << "\n";
                break;
            case DataType::Tag::Character:
                if (field.rank == 0)
                    w << "    " << var_f << " = TRANSFER(" << var_c << ", " << var_f << " )\n";
                break;
            case DataType::Tag::Derived:
                break;
            }
        }
    }

    w << " END SUBROUTINE " << routine_name << "\n\n";
}

void gen_copy_f2c(std::ostream &w, const Module &mod, const DataType::Derived &ddt)
{
    std::string routine_name = mod.nickname + "_F2C_Copy" + ddt.name_short;

    w << " SUBROUTINE " << routine_name << "( " << ddt.name_short
      << "Data, ErrStat, ErrMsg, SkipPointers  )\n";
    w << "    TYPE(" << ddt.type_fortran << "), INTENT(INOUT) :: " << ddt.name_short << "Data\n";
    w << "    INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n";
    w << "    CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n";
    w << "    LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers\n";
    w << "    ! \n";
    w << "    LOGICAL                        :: SkipPointers_local\n";
    w << "    ErrStat = ErrID_None\n";
    w << "    ErrMsg  = \"\"\n\n";
    w << "    IF (PRESENT(SkipPointers)) THEN\n";
    w << "       SkipPointers_local = SkipPointers\n";
    w << "    ELSE\n";
    w << "       SkipPointers_local = .false.\n";
    w << "    END IF\n";

    for (const auto &field : ddt.fields)
    {

        // If field doesn't have a data type, continue
        if (field.data_type == nullptr)
        {
            continue;
        }

        // If field is a derived type, print warning and continue
        if (field.data_type->tag == DataType::Tag::Derived)
        {
            std::cerr << "Registry WARNING: derived data type " << field.name << " of type "
                      << field.data_type->derived.name << " is not passed through F-C interface\n";
            continue;
        }

        std::string var_f = ddt.name_short + "Data%" + field.name;
        std::string var_c = ddt.name_short + "Data%C_obj%" + field.name;

        if (field.is_pointer)
        {
            w << "\n    ! -- " << field.name << " " << ddt.name_short << " Data fields\n";
            w << "    IF ( .NOT. SkipPointers_local ) THEN\n";
            w << "       IF ( .NOT. ASSOCIATED(" << var_f << ")) THEN \n";
            w << "          " << var_c << "_Len = 0\n";
            w << "          " << var_c << " = C_NULL_PTR\n";
            w << "       ELSE\n";
            w << "          " << var_c << "_Len = SIZE(" << var_f << ")\n";
            w << "          IF (" << var_c << "_Len > 0) &\n";
            w << "             " << var_c << " = C_LOC( " << var_f << "(";
            for (int d = 1; d <= field.rank; d++)
            {
                w << (d > 1 ? "," : "") << " LBOUND(" << var_f << "," << d << ")";
            }
            w << " ) )\n";
            w << "       END IF\n";
            w << "    END IF\n";
        }
        else if (!field.is_allocatable)
        {

            switch (field.data_type->tag)
            {
            case DataType::Tag::Real:
            case DataType::Tag::Integer:
            case DataType::Tag::Logical:
                w << "    " << var_c << " = " << var_f << "\n";
                break;
            case DataType::Tag::Character:
                if (field.rank == 0)
                    w << "    " << var_c << " = TRANSFER(" << var_f << ", " << var_c << " )\n";
                break;
            case DataType::Tag::Derived:
                break;
            }
        }
    }

    w << " END SUBROUTINE " << routine_name << "\n\n";
}