#include <fstream>

#include "registry.hpp"
#include "templates.hpp"

void output_template(std::string &module_name, std::string &module_nickname, bool overwrite,
                     bool is_template);

const std::string usage_template = R""""(
Usage: openfast_registry registryfile [options] -or-
          [-force] [-template|-registry] ModuleName ModName 
Options:
    -h                this summary
    -I <dir>          look for usefrom files in directory "dir"
    -O <dir>          generate types files in directory "dir"
    -incsubs          generate the pack/unpack/copy/destroy subroutines to be included in another file
    -noextrap         do not generate ModName_Input_ExtrapInterp or ModName_Output_ExtrapInterp routines
    -D<SYM>           define symbol for conditional evaluation inside registry file
    -ccode            generate additional code for interfacing with C/C++
    -keep             do not delete temporary files from registry program
    -shownodes        output a listing of the nodes in registry's AST
  === alternate usage for generating templates ===
    -template ModuleName ModName
                 Generate a template Module file none exists
    -registry ModuleName ModName
                 Generate a template registry file if none exists
    -force Force generating of template or registry file
  (the / character can be used in place of - when specifying options)
)"""";

int main(int argc, char *argv[])
{
    std::cerr << std::endl;
    std::cerr << "------------------------------------------------------------" << std::endl;
    std::cerr << "-------------------- OpenFAST Registry ---------------------" << std::endl;
    std::cerr << "------------------------------------------------------------" << std::endl;

    // Read command line arguments into a vector
    std::vector<std::string> arguments;
    for (int i = 0; i < argc; ++i)
    {
        arguments.push_back(argv[i]);
    }

    std::string out_dir = "."; // if no OutDir is listed, use current directory
    std::string inp_file_path;
    std::string module_name, module_nickname;
    bool output_force_template = false;

    // Create registry object
    Registry reg;

    // Loop through arguments
    for (auto it = arguments.begin(); it != arguments.end(); ++it)
    {
        auto arg = *it;

        if ((arg.compare("-force") == 0) || (arg.compare("/force") == 0))
        {
            output_force_template = true;
        }
        else if ((arg.compare("-ccode")) == 0 || (arg.compare("/ccode")) == 0)
        {
            reg.gen_c_code = true;
        }
        else if ((arg.compare("-noextrap")) == 0 || (arg.compare("/noextrap")) == 0)
        {
            reg.no_extrap_interp = true;
        }
        else if ((arg.compare("-shownodes")) == 0 || (arg.compare("/shownodes")) == 0)
        {
        }
        else if ((arg.compare("-O")) == 0 || (arg.compare("/O")) == 0)
        {
            std::advance(it, 1);
            if (it != arguments.end())
            {
                out_dir = *it;
            }
        }
        else if ((arg.compare("-I")) == 0 || (arg.compare("/I")) == 0)
        {
            std::advance(it, 1);
            if (it != arguments.end())
            {
                reg.include_dirs.push_back(*it);
            }
        }
        else if ((arg.compare("-incsubs")) == 0 || (arg.compare("/incsubs")) == 0)
        {
            reg.gen_inc_subs = true;
        }
        else if ((arg.compare("-template")) == 0 || (arg.compare("-registry")) == 0 ||
                 (arg.compare("/template")) == 0 || (arg.compare("/registry")) == 0)
        {
            std::advance(it, 1);
            if (it != arguments.end())
            {
                module_name = *it;
            }
            else
            {
                std::cerr << usage_template;
                return EXIT_FAILURE;
            }
            std::advance(it, 1);
            if (it != arguments.end())
            {
                module_nickname = *it;
            }
            else
            {
                std::cerr << usage_template;
                return EXIT_FAILURE;
            }

            bool is_template = arg.substr(1).compare("template") == 0;

            output_template(module_name, module_nickname, output_force_template, is_template);
            return EXIT_SUCCESS;
        }
        else if ((arg.compare("-h") == 0) || (arg.compare("/h") == 0))
        {
            std::cerr << usage_template;
            return EXIT_SUCCESS;
        }
        else
        {
            // Set input file path
            inp_file_path = arg;

            // Replace backslashes with forward slashes in path
            std::string path = std::regex_replace(arg, std::regex("\\\\"), "/");

            // If path contains / remove everything after it
            auto slash_index = path.find_last_of("/");
            if (slash_index != std::string::npos)
                path = path.substr(0, slash_index);

            // Add input file directory to list of include directories
            reg.include_dirs.push_back(path);
        }
    }

    // If input file name was not specified, exit with error
    if (inp_file_path.empty())
    {
        std::cerr << usage_template;
        return EXIT_FAILURE;
    }

    // Parse the registry file
    reg.parse(inp_file_path, 0);

    // Generate module files
    reg.gen_module_files(out_dir);
}

void output_template(std::string &module_name, std::string &module_nickname, bool overwrite,
                     bool is_template)
{
    // Create file name depending on if template or registry
    std::string fname = module_name + (is_template ? ".f90" : "_Registry.txt");

    // If overwrite not requested and file exists, return error
    if (!overwrite)
    {
        std::ifstream infile(fname);
        if (infile.good())
        {
            std::cerr << "Registry exiting. Attempt to overwrite file (" << fname;
            std::cerr << ") . Move out of the way or specify -force before -template option. "
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    // Open output file, return on error
    std::ofstream outfile(fname);
    if (!outfile.is_open())
    {
        std::cerr << "Registry exiting. Failure opening " << fname << std::endl;
        exit(EXIT_FAILURE);
    }

    // Select file contents
    auto contents = (is_template ? module_template : registry_template);

    // Populate module name and module nickname
    contents = std::regex_replace(contents, std::regex("ModuleName"), module_name);
    contents = std::regex_replace(contents, std::regex("ModName"), module_nickname);

    // Output contents to file
    outfile << contents;

    std::cerr << "Created " << (is_template ? "template" : "registry") << " file '" << fname << "'" << std::endl;
}
