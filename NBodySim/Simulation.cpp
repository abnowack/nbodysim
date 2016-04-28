#include "Simulation.h"
#include <fstream>
#include <string>

std::array<double, 3> Simulation::gravity_accel(std::array<double, 3> r, std::array<double, 3> r_other, double mass_other) {
	//double G = -2.95912208286e-4;
	double G = -39.48763963;
	double r_norm = std::sqrt(std::pow(r[0] - r_other[0], 2.) +
		std::pow(r[1] - r_other[1], 2.) +
		std::pow(r[2] - r_other[2], 2.));

	decltype(r) accel;
	accel[0] = G * mass_other * (r[0] - r_other[0]) / std::pow(r_norm, 3.);
	accel[1] = G * mass_other * (r[1] - r_other[1]) / std::pow(r_norm, 3.);
	accel[2] = G * mass_other * (r[2] - r_other[2]) / std::pow(r_norm, 3.);

	return accel;
}

void Simulation::add_particle(rv3 &rv, double mass, std::string name, bool interacting) {
	rv3 dummy; // ugh, I dont know TODO: remove dummy
	dummy.rv[0] = 0;
	dummy.rv[1] = 0;
	dummy.rv[2] = 0;
	dummy.rv[3] = 0;
	dummy.rv[4] = 0;
	dummy.rv[5] = 0;
	rv_.push_back(rv);
	rv_new_.emplace_back(dummy);
	mass_.push_back(mass);
	name_.emplace_back(name);
	id_.push_back(id_count_);
	if (interacting) {
		unsigned int index = static_cast<unsigned int>(rv_.size() - 1);
		interacting_index_.push_back(index);
	}

	id_count_++;
}

std::array<double, 3> Simulation::calc_accel(unsigned int i, bool use_new) {
	rv3 &r = use_new ? rv_new_[i] : rv_[i];

	std::array<double, 3> total_accel = { 0., 0., 0. };
	for (unsigned int index_i = 0; index_i < interacting_index_.size(); index_i++) {
		auto index = interacting_index_[index_i];
		if (id_[i] == id_[index]) continue;
		auto j_accel = gravity_accel(r.r, rv_[index].r, mass_[index]);
		total_accel[0] += j_accel[0];
		total_accel[1] += j_accel[1];
		total_accel[2] += j_accel[2];
	}
	return total_accel;
}

void Simulation::apply_pefrl(unsigned int i) {
	double xi = +0.1786178958448091E+00;
	double lambda = -0.2123418310626054E+00;
	double chi = -0.6626458266981849E-01;

	rv3 &rv_new = rv_new_[i];
	rv3 &rv = rv_[i];

	// I
	for (unsigned int j = 0; j < 3; j++) {
		rv_new.r[j] = rv.r[j] + xi * step_ * rv.v[j];
	}
	auto a = calc_accel(i, true);
	for (unsigned int j = 0; j < 3; j++) {
		rv_new.v[j] = rv.v[j] + (1 - 2. * lambda) * step_ / 2. * a[j];
	}
	// II
	for (unsigned int j = 0; j < 3; j++) {
		rv_new.r[j] += chi * step_ * rv_new.v[j];
	}
	a = calc_accel(i, true);
	for (unsigned int j = 0; j < 3; j++) {
		rv_new.v[j] += lambda * step_ * a[j];
	}
	// III
	for (unsigned int j = 0; j < 3; j++) {
		rv_new.r[j] += (1. - 2. * (chi + xi)) * step_ * rv_new.v[j];
	}
	a = calc_accel(i, true);
	for (unsigned int j = 0; j < 3; j++) {
		rv_new.v[j] += lambda * step_ * a[j];
	}
	// IV
	for (unsigned int j = 0; j < 3; j++) {
		rv_new.r[j] += chi * step_ * rv_new.v[j];
	}
	// v(t + h)
	a = calc_accel(i, true);
	for (unsigned int j = 0; j < 3; j++) {
		rv_new.v[j] += (1 - 2 * lambda) * step_ / 2. * a[j];
	}
	// r(t + h)
	for (unsigned int j = 0; j < 3; j++) {
		rv_new.r[j] += xi * step_ * rv_new.v[j];
	}
}


void Simulation::update() {
	for (unsigned int i = 0; i < rv_.size(); i++) {
		apply_pefrl(i);
	}
	for (unsigned int i = 0; i < rv_.size(); i++) {
		rv_[i] = rv_new_[i];
	}

	time_ += step_;
}

void Simulation::remove_global_momentum() {
	std::array<double, 3> total_mom = { 0, 0, 0 };
	double total_mass = 0;
	for (unsigned int i = 0; i < rv_.size(); i++) {
		total_mom[0] += mass_[i] * rv_[i].v[0];
		total_mom[1] += mass_[i] * rv_[i].v[1];
		total_mom[2] += mass_[i] * rv_[i].v[2];
		total_mass += mass_[i];
	}
	// remove total_mom / size(particles_) momentum from each particle
	for (unsigned int i = 0; i < rv_.size(); i++) {
		rv_[i].v[0] -= total_mom[0] / total_mass;
		rv_[i].v[1] -= total_mom[1] / total_mass;
		rv_[i].v[2] -= total_mom[2] / total_mass;
	}
}


void Simulation::add_particles_csv(std::string fname) {
	std::string name;
	rv3 rv;
	double mass;

	std::ifstream csv_file;
	std::string read_data;
	csv_file.open(fname);

	while (csv_file.good()) {
		getline(csv_file, name, ',');

		getline(csv_file, read_data, ',');
		mass = std::stod(read_data);

		getline(csv_file, read_data, ',');
		rv.r[0] = std::stod(read_data);

		getline(csv_file, read_data, ',');
		rv.r[1] = std::stod(read_data);

		getline(csv_file, read_data, ',');
		rv.r[2] = std::stod(read_data);

		getline(csv_file, read_data, ',');
		rv.v[0] = std::stod(read_data);

		getline(csv_file, read_data, ',');
		rv.v[1] = std::stod(read_data);

		getline(csv_file, read_data);
		rv.v[2] = std::stod(read_data);

		if (csv_file.good() && (mass > 0)) {
			add_particle(rv, mass, name, true);
		}
		else if (csv_file.good()) {
			add_particle(rv, mass, name, false);
		}
	}

	csv_file.close();
}
