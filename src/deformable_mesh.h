#pragma once

#include <string>
#include <Eigen/Core>

namespace mti855 {
namespace physics {

	enum class deformation_type_t
	{
		rotation,
		linear,
		quadratic
	};

	/* https://www.cs.drexel.edu/~david/Classes/Papers/MeshlessDeformations_SIG05.pdf */
	class deformable_mesh
	{
	public:
		Eigen::MatrixX3d const& V() const { return x_; }
		Eigen::MatrixX3i const& F() const { return F_; }

		Eigen::MatrixX3d const& forces() const { return f_; }
		Eigen::MatrixX3d const& velocities() const { return v_; }
		Eigen::MatrixX3d const& positions() const { return x_; }

		Eigen::MatrixX3d& V() { return x_; }
		Eigen::MatrixX3i& F() { return F_; }

		bool load_from_obj(std::string const& filepath);
		
		void apply_gravity() { f_.rowwise() += Eigen::RowVector3d{ 0.0, -9.8, 0.0 }; }
		Eigen::MatrixX3d apply_force(Eigen::Vector3d const& center, Eigen::Vector3d const& force);
		void integrate(double dt);
		void integrate_quadratic(double dt);

		void set_rayleigh_beta(double Rb) { Rb_ = Rb; }
		void set_rayleigh_alpha(double Ra) { Ra_ = Ra; }

		double rayleigh_beta() const { return Rb_; }
		double rayleigh_alpha() const { return Ra_; }

		void set_beta(double beta) { beta_ = beta; }
		double beta() const { return beta_; }

		double tau() const { return tau_; }
		void set_tau(double tau) { tau_ = tau; }

		double perturbation() const { return perturbation_;  }
		void set_perturbation(double perturbation) { perturbation_ = perturbation; }

		deformation_type_t deformation_type() const { return deformation_type_; }
		void set_deformation_type(deformation_type_t type) { deformation_type_ = type; }

		bool is_fixed(unsigned int vi) const { return fixed_(vi); }
		void set_fixed(unsigned int vi, bool fixed = true) { fixed_(vi) = fixed; }
		unsigned int count_fixed() const { return fixed_.count(); }
		void clear_fixed() { fixed_.fill(false); }

		std::size_t get_memory_in_bytes() const
		{
			std::size_t mem = 0;
			mem += V0_.rows() * V0_.cols() * sizeof(decltype(V0_)::CoeffReturnType);
			mem += Q_.rows()  * Q_.cols()  * sizeof(decltype(Q_)::CoeffReturnType);
			mem += AqqInv_.rows() * AqqInv_.cols() * sizeof(decltype(AqqInv_)::CoeffReturnType);
			mem += Qquadratic_.rows() * Qquadratic_.cols() * sizeof(decltype(Qquadratic_)::CoeffReturnType);
			mem += AqqInvQuadratic_.rows() * AqqInvQuadratic_.cols() * sizeof(decltype(AqqInvQuadratic_)::CoeffReturnType);
			mem += x_.rows() * x_.cols() * sizeof(decltype(x_)::CoeffReturnType);
			mem += v_.rows() * v_.cols() * sizeof(decltype(v_)::CoeffReturnType);
			mem += fixed_.rows() * fixed_.cols() * sizeof(decltype(fixed_)::CoeffReturnType);
			mem += f_.rows() * f_.cols() * sizeof(decltype(f_)::CoeffReturnType);
			mem += F_.rows() * F_.cols() * sizeof(decltype(F_)::CoeffReturnType);

			return mem;
		}

	private:
		void reset();

		// rest space vertices
		Eigen::MatrixX3d V0_;

		// matrix of qi = xi0 - xicm0 for linear deformations
		Eigen::MatrixX3d Q_;

		// precomputed for linear deformations
		Eigen::Matrix3d AqqInv_;

		// matrix of qi for quadratic deformations
		Eigen::Matrix<double, Eigen::Dynamic, 9> Qquadratic_;

		// precomputed for quadratic deformations
		Eigen::Matrix<double, 9, 9> AqqInvQuadratic_;

		/* state vectors */
		// positions
		Eigen::MatrixX3d x_;
		// velocities
		Eigen::MatrixX3d v_;

		Eigen::Array<bool, Eigen::Dynamic, 1> fixed_;

		/* forces */
		Eigen::MatrixX3d f_;

		// triangles
		Eigen::MatrixX3i F_;

		// Rayleigh damping
		double Rb_ = 0.0;
		double Ra_ = 0.0;

		// elasticity params
		double tau_ = 1.0;

		// Beta
		double beta_ = 0.0;

		// regularization perturbation
		double perturbation_ = 0.1;

		deformation_type_t deformation_type_;
	};

} // physics
} // mti855