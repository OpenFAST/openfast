//********************************************************************************************************************************
// LICENSING
// Copyright (C) 2023-2025  National Renewable Energy Laboratory
// Copyright (C) 2025-2026  National Laboratory of the Rockies
//
//    This file is part of the OpenFAST Registry.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//********************************************************************************************************************************
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