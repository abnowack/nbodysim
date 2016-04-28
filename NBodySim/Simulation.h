#include <array>
#include <vector>
#include <string>

union rv3 {
	std::array<double, 6> rv;
	struct {
		std::array<double, 3> r;
		std::array<double, 3> v;
	};
};

class Simulation {
public:
	Simulation(double step) : step_(step)
	{
	}

	static std::array<double, 3> gravity_accel(std::array<double, 3> r, std::array<double, 3> r_other, double mass_other);

	void add_particle(rv3 &rv, double mass, std::string name, bool interacting = false);
	std::array<double, 3> calc_accel(unsigned int i, bool use_new = false);
	void apply_pefrl(unsigned int i);
	void update();
	void remove_global_momentum();
	void add_particles_csv(std::string fname);

	std::vector<rv3> rv_;
	std::vector<rv3> rv_new_;
	std::vector<double> mass_;
	std::vector<std::string> name_;
	std::vector<unsigned int> id_;
	std::vector<unsigned int> interacting_index_;

private:
	double step_ = 1;
	double time_ = 0;
	unsigned int id_count_ = 0;
};

