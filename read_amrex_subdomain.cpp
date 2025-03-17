
#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace
{
    struct IntVect
    {
	int a[3];
    };

    std::ostream& operator<< (std::ostream& os, IntVect const& iv)
    {
	os << '(' << iv.a[0] << ',' << iv.a[1] << ',' << iv.a[2] << ')';
	return os;
    }

    std::istream& operator>> (std::istream& is, IntVect& iv)
    {
	is.ignore(1000000, '(');
	is >> iv.a[0];
	is.ignore(1000000, ',');
	is >> iv.a[1];
	is.ignore(1000000, ',');
	is >> iv.a[2];
	is.ignore(1000000, ')');
	return is;
    }

    struct Box
    {
	IntVect lo, hi;
    };

    std::ostream& operator<< (std::ostream& os, Box const& box)
    {
	os << '(' << box.lo << ' ' << box.hi << ')';
	return os;
    }

    std::istream& operator>> (std::istream& is, Box& box)
    {
	is.ignore(1000000, '(');
	is >> box.lo >> box.hi;
	is.ignore(1000000, ')');
	is.ignore(1000000, ')');
	return is;
    }

    struct FOD
    {
	std::string file;
	std::size_t offset;
    };

    std::vector<Box> grids;
    std::vector<FOD> fods;

    void read_file (std::string const& File, std::vector<char>& file_char_vec)
    {
	std::ifstream ifs(File, std::ios::in);
	if (!ifs.good()) {
	    std::cout << "Failed to read " << File << "\n";
	    std::abort();
	}
	ifs.seekg(0, std::ios::end);
	auto file_length = static_cast<std::streamoff>(ifs.tellg());
	ifs.seekg(0, std::ios::beg);
	file_char_vec.resize(file_length+1);
	ifs.read(file_char_vec.data(), file_length);
	file_char_vec[file_length] = '\0';
    }
}

extern "C"
{

void read_amrex_header (char const* name, int* dims, double* origin, double* dx,
			double* time)
{
    std::vector<char> file_char_vec;
    read_file(std::string(name)+"/Header", file_char_vec);

    std::array<double,3> prob_lo;
    {
	std::istringstream is(std::string(file_char_vec.data()), std::istringstream::in);

	std::string file_version;
	is >> file_version;

	int ncomp;
	is >> ncomp;
	assert(ncomp == 3);

	is.ignore(1000000, '\n');
	for (int i = 0; i < ncomp; ++i) {
	    std::string tmp;
	    std::getline(is, tmp); // variable names;
	}

	int spacedim, finest_level;
	is >> spacedim >> *time >> finest_level;
	assert(spacedim == 3 && finest_level == 0);

	int nlevels = finest_level+1;

	std::array<double,3> prob_hi;

	for (int i = 0; i < spacedim; ++i) {
	    is >> prob_lo[i];
	}
	for (int i = 0; i < spacedim; ++i) {
	    is >> prob_hi[i];
	}

	is.ignore(1000000, '\n');
    
	for (int i = 0; i < nlevels; ++i) {
	    std::string tmp;
	    std::getline(is, tmp); // domain box
	    is.ignore(1000000, '\n');
	}

	int level_steps;
	for (int i = 0; i < nlevels; ++i) {
	    is >> level_steps;
	}
	std::cout << " level_steps = " << level_steps << std::endl;

	for (int ilev = 0; ilev < nlevels; ++ilev) {
	    is >> dx[0] >> dx[1] >> dx[2];
	}
    }

    read_file(std::string(name)+"/Level_0/Cell_H", file_char_vec);

    {
	std::istringstream is(std::string(file_char_vec.data()), std::istringstream::in);
	int version, how, ncomp, ng;
	is >> version >> how >> ncomp >> ng;
	assert(ncomp == 3 && ng == 0);

	// Read BoxArray
	int nboxes, tmp;
	is.ignore(1000000, '(') >> nboxes >> tmp;
	grids.resize(nboxes);
	for (auto& b : grids) {
	    is >> b;
	}
	is.ignore(1000000, ')');

	int lo[3] = {std::numeric_limits<int>::max(),
	             std::numeric_limits<int>::max(),
	             std::numeric_limits<int>::max()};
	int hi[3] = {std::numeric_limits<int>::lowest(),
	             std::numeric_limits<int>::lowest(),
	             std::numeric_limits<int>::lowest()};
	for (auto const& b : grids) {
	    for (int i = 0; i < 3; ++i) {
		lo[i] = std::min(lo[i], b.lo.a[i]);
		hi[i] = std::max(hi[i], b.hi.a[i]);
	    }
	}
	for (int i = 0; i < 3; ++i) {
	    dims[i] = hi[i] - lo[i] + 1;
	    origin[i] = prob_lo[i] + i*dx[i];
	}

	// Read FabOnDisk
	fods.resize(nboxes);
    }
}

}
