#include <iostream>

#include "H5Cpp.h"
#include "Simulation.h"

#pragma warning(disable:4996) // security flags

// TODO: Write results to HDF5 file
// TODO: Have Simulation write HDF5 automatically
// TODO: Read particle information via json file on commandline
// TODO: Create Python animated gif maker of particle interactions

class SimulationWriter {
public:

	SimulationWriter(Simulation &sim, std::string fname = "out.h5", unsigned int steps = 10) {
		auto nparticles = sim.rv_.size();

		H5::H5File file(fname.c_str(), H5F_ACC_TRUNC);
		hsize_t fdims[3] = { steps, nparticles, 6 };
		H5::DataSpace pos_dataspace = H5::DataSpace(3, fdims);

		xyz_dataset_ = file.createDataSet("data", H5::PredType::NATIVE_DOUBLE, pos_dataspace);
		xyz_dataspace_ = xyz_dataset_.getSpace();

		// Write Particles Names
		char cname[256];
		hsize_t adims[1] = { nparticles };
		H5::DataSpace info_d(1, adims);

		H5::StrType strdatatype(H5::PredType::C_S1, 256);
		auto info_dataset = file.createDataSet("name", strdatatype, info_d);
		auto info_dataspace = info_dataset.getSpace();
		hsize_t aoffset[1] = { 1 }, acount[1] = { 1 }, astride[1] = { 1 }, ablock[1] = { 1 }, adimsm[1] = { 1 };
		auto info_mspace = H5::DataSpace(1, adimsm, NULL);
		for (unsigned int i = 0; i < nparticles; i++) {
			aoffset[0] = i;
			info_dataspace.selectHyperslab(H5S_SELECT_SET, acount, aoffset, astride, ablock);
			strcpy_s(cname, sim.name_[i].c_str());
			H5std_string strwritebuf = cname;
			info_dataset.write((sim.name_[i] + "\0").c_str(), strdatatype, info_mspace, info_dataspace);
		}

		// Write Particle Masses
		hsize_t mdims[1] = { nparticles };
		H5::DataSpace mass_d(1, mdims);
		auto mass_dataset = file.createDataSet("mass", H5::PredType::NATIVE_DOUBLE, mass_d);
		auto mass_dataspace = mass_dataset.getSpace();
		hsize_t moffset[1] = { 1 }, mcount[1] = { 1 }, mstride[1] = { 1 }, mblock[1] = { 1 }, mdimsm[1] = { 1 };
		auto mass_mspace = H5::DataSpace(1, mdimsm, NULL);
		for (unsigned int i = 0; i < nparticles; i++) {
			moffset[0] = i;
			mass_dataspace.selectHyperslab(H5S_SELECT_SET, mcount, moffset, mstride, mblock);
			mass_dataset.write(&sim.mass_[i], H5::PredType::NATIVE_DOUBLE, mass_mspace, mass_dataspace);
		}

		count_[1] = nparticles;
		dimsm_[0] = nparticles;
		memspace_ = H5::DataSpace(2, dimsm_, NULL);
	}

	void writeParticleData(unsigned int istep, Simulation &sim) {
		offset_[0] = istep;
		xyz_dataspace_.selectHyperslab(H5S_SELECT_SET, count_, offset_, stride_, block_);
		xyz_dataset_.write(&sim.rv_[0], H5::PredType::NATIVE_DOUBLE, memspace_, xyz_dataspace_);
	}

	// HDF5 stuff
	H5::DataSet xyz_dataset_;
	hsize_t dimsm_[2] = { 1, 6 }, offset_[3] = { 0, 0, 0 }, count_[3] = { 1, 1, 6 }, stride_[3] = { 1, 1, 1 }, block_[3] = { 1, 1, 1 };
	H5::DataSpace memspace_, xyz_dataspace_;;
};

int main(int argc, char* argv[])
{
	Simulation sim(0.05);
	//sim.add_particles_csv("H:\\nbody\\oldobjects.csv");
	sim.add_particles_csv("H:\\nbody\\objects.csv");
	std::cout << sim.rv_.size() << " " << sim.interacting_index_.size() << std::endl;
	sim.remove_global_momentum();
	unsigned int nsteps = 2000000, write_nsteps = 200;
	unsigned int total_write_steps = nsteps / write_nsteps;
	SimulationWriter simwrite(sim, "perturbed.h5", total_write_steps);

	unsigned int write_steps = 0;
	for (unsigned int i = 0; i < nsteps; i++) {
		if ((i % write_nsteps) == 0) {
			simwrite.writeParticleData(write_steps, sim);
			write_steps++;
		}
		sim.update();
		if (i % (nsteps / 100) == 0) {
			std::cout << i << " / " << nsteps << std::endl;
		}
	}

	simwrite.xyz_dataset_.close();

	std::cout << "done" << std::endl;
	std::cin.get();

	return 0;
}
