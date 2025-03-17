
#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
    int coarsen (int i, int ratio)
    {
	return (i<0) ? -std::abs(i+1)/ratio-1 : i/ratio;
    }

    struct IntVect
    {
	int a[3];
	int& operator[] (int i) { return a[i]; }
	int const& operator[] (int i) const { return a[i]; }

	bool operator== (IntVect const& rhs) const noexcept
	{
	    return (a[0] == rhs.a[0])
		&& (a[1] == rhs.a[1])
		&& (a[2] == rhs.a[2]);
	}

	struct shift_hasher
	{
	    std::size_t operator() (IntVect const& vec) const noexcept
	    {
                static constexpr unsigned shift1 = sizeof(size_t)>=8 ? 20 : 10;
                static constexpr unsigned shift2 = sizeof(size_t)>=8 ? 40 : 20;
                return static_cast<std::size_t>(vec[0]) ^
                      (static_cast<std::size_t>(vec[1]) << shift1) ^
                      (static_cast<std::size_t>(vec[2]) << shift2);
            }
        };
    };

    std::ostream& operator<< (std::ostream& os, IntVect const& iv)
    {
	os << '(' << iv[0] << ',' << iv[1] << ',' << iv[2] << ')';
	return os;
    }

    std::istream& operator>> (std::istream& is, IntVect& iv)
    {
	is.ignore(1000000, '(');
	is >> iv[0];
	is.ignore(1000000, ',');
	is >> iv[1];
	is.ignore(1000000, ',');
	is >> iv[2];
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

    struct HashMap
    {
	IntVect blo;
	IntVect bhi;
	IntVect maxext;
	std::unordered_map<IntVect,std::vector<int>,IntVect::shift_hasher> hash;
    };

    std::vector<Box> grids;
    std::vector<FOD> fods;
    HashMap hashmap;

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

	for (int i = 0; i < 3; ++i) {
	    hashmap.blo[i] = std::numeric_limits<int>::max();
	    hashmap.bhi[i] = std::numeric_limits<int>::lowest();
	    hashmap.maxext[i] = 1;
	}
	for (auto const& b : grids) {
	    for (int i = 0; i < 3; ++i) {
		hashmap.blo[i] = std::min(hashmap.blo[i], b.lo[i]);
		hashmap.bhi[i] = std::max(hashmap.bhi[i], b.hi[i]);
		hashmap.maxext[i] = std::max(hashmap.maxext[i], b.hi[i]-b.lo[i]+1);
	    }
	}
	for (int i = 0; i < 3; ++i) {
	    dims[i] = hashmap.bhi[i] - hashmap.blo[i] + 1;
	    origin[i] = prob_lo[i] + i*dx[i];
	}

	// Read FabOnDisk
	int nfabs;
	is >> nfabs;
	assert(nfabs == nboxes);
	fods.resize(nfabs);
	for (int i = 0; i < nfabs; ++i) {
	    std::string stmp, str;
	    is >> stmp >> str >> fods[i].offset;
	    fods[i].file = std::string(name) + "/Level_0/" + str;
	}
    }

    // Build hash
    auto nboxes = int(grids.size());
    for (int i = 0; i < nboxes; ++i) {
	IntVect key;
	for (int idim = 0; idim < 3; ++idim) {
	    key[idim] = coarsen(grids[i].lo[idim], hashmap.maxext[idim]);
	}
	hashmap.hash[key].push_back(i);
    }
}

void read_amrex_subdomain (double* a, int const* a_lo, int const* a_hi)
{
    long long istride_a = 3;
    long long jstride_a = istride_a * (a_hi[0]-a_lo[0]+1);
    long long kstride_a = jstride_a * (a_hi[1]-a_lo[1]+1);

    std::map<std::string,std::ifstream> ifs_map;

    IntVect clo, chi, alo, ahi;
    for (int idim = 0; idim < 3; ++idim) {
	alo[idim] = a_lo[idim] + hashmap.blo[idim];
	ahi[idim] = a_hi[idim] + hashmap.blo[idim]; // Yes blo is the offset.
	clo[idim] = coarsen(std::max(alo[idim],hashmap.blo[idim]),
			    hashmap.maxext[idim]);
	chi[idim] = coarsen(std::min(ahi[idim],hashmap.bhi[idim]),
			    hashmap.maxext[idim]);
    }

    // -1 because we use the lower corner for hash
    for (int kk = clo[2]-1; kk <= chi[2]; ++kk) {
    for (int jj = clo[1]-1; jj <= chi[1]; ++jj) {
    for (int ii = clo[0]-1; ii <= chi[0]; ++ii) {
	IntVect key{ii,jj,kk};
	auto it = hashmap.hash.find(key);
	if (it != hashmap.hash.cend()) {
	    for (int ibox : it->second) {
		int ilo[3], ihi[3];
		bool ok = true;
		for (int idim = 0; idim < 3; ++idim) {
		    ilo[idim] = std::max(alo[idim], grids[ibox].lo[idim]);
		    ihi[idim] = std::min(ahi[idim], grids[ibox].hi[idim]);
		    ok = ok && (ilo[idim] <= ihi[idim]);
		}

		if (ok) {
		    auto& ifs = ifs_map[fods[ibox].file];
		    if (! ifs.is_open()) {
			ifs.open(fods[ibox].file, std::ios::in|std::ios::binary);
		    }
		    ifs.seekg(fods[ibox].offset, std::ios::beg);
		    char c;
		    bool badfab = false;
		    ifs >> c;
		    if(c != 'F') {
			badfab = true;
		    }
		    ifs >> c;
		    if(c != 'A') {
			badfab = true;
		    }
		    ifs >> c;
		    if (c != 'B') {
			badfab = true;
		    }
		    if (badfab) {
			std::cout << "Bad Fab in " << fods[ibox].file << "\n";
			std::abort();
		    }
		    for (int i = 0; i < 5; ++i) { // Real descriptor has 5 ')'s
			ifs.ignore(1000000, ')');
		    }
		    Box box;
		    ifs >> box;
		    bool badbox = false;
		    for (int idim = 0; idim < 3; ++idim) {
			if ((box.lo[idim] != grids[ibox].lo[idim]) ||
			    (box.hi[idim] != grids[ibox].hi[idim])) {
			    badbox = true;
			}
		    }
		    if (badbox) {
			std::cout << "Bad box in " << fods[ibox].file << ", " << box
				  << ", " << grids[ibox] << "\n";
			std::abort();
		    }
		    int ncomp;
		    ifs >> ncomp;
		    assert(ncomp == 3);
		    ifs.ignore(1000000, '\n');

		    long long jstride_p = box.hi[0]-box.lo[0]+1;
		    long long kstride_p = jstride_p * (box.hi[1]-box.lo[1]+1);
		    long long nstride_p = kstride_p * (box.hi[2]-box.lo[2]+1);
		    auto nreals = nstride_p * ncomp;
		    auto nbytes = sizeof(double)*nreals;
		    auto* p = (double*)std::malloc(nbytes);
		    ifs.read((char*)p, std::streamsize(nbytes));

		    for (int k = ilo[2]; k <= ihi[2]; ++k) {
			int ka = k - alo[2];
			int kp = k - box.lo[2];
			for (int j = ilo[1]; j <= ihi[1]; ++j) {
			    int ja = j - alo[1];
			    int jp = j - box.lo[1];
			    for (int i = ilo[0]; i <= ihi[0]; ++i) {
				int ia = i - alo[0];
				int ip = i - box.lo[0];
				long long aoff = ia*istride_a + ja*jstride_a + ka*kstride_a;
				long long poff = ip + jp*jstride_p + kp*kstride_p;
				a[aoff  ] = p[poff            ];
				a[aoff+1] = p[poff+nstride_p  ];
				a[aoff+2] = p[poff+nstride_p*2];
			    }
			}
		    }

		    std::free(p);
		}
	    }
	}}}
    }
}

}
