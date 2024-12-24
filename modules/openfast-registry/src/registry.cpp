#include <cstdlib>
#include <cctype>

#include "registry.hpp"

void Registry::gen_module_files(std::string const &out_dir)
{
    // Find root module
    std::shared_ptr<Module> mod;
    for (auto &it : this->modules)
    {
        if (it.second->is_root)
        {
            mod = it.second;
            break;
        }
    }

    // If module not found, return error
    if (mod == nullptr)
    {
        std::cerr << "unable to find root module" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Generate fortran module
    this->gen_fortran_module(*mod, out_dir);

    // Generate C code
    if (this->gen_c_code)
        this->gen_c_module(*mod, out_dir);
}

std::string tolower(std::string s)
{
    for (auto &c : s)
        c = std::tolower(c);
    return s;
}