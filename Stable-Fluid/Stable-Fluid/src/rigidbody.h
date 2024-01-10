#ifndef  RIGIDBODY_H
#define RIGIDBODY_H

#ifndef PI
#define PI 3.1415926f
#endif
#ifndef G_
#define G_ 9.8f * 5
#endif
#ifndef Gamma_
#define Gamma_ 0.6f
#endif

#include "field.h"
#include "geometry.h"
#include "util.h"

#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class RigidBody {
private:
	// const
	std::vector<glm::vec3> vertices_;
	std::vector<Face> faces_;
	float M_;
	glm::mat3 I_ref_;
	glm::mat3 I_inv_ref_;
	glm::vec3 c0_;
	// var
	glm::vec3 c_;
	glm::quat rotationQuaternion_;
	glm::vec3 Vc_;
	glm::vec3 omega_;
	glm::vec3 L_;// L_ = I() * omega_
	// collision
	float e_ = 0.3f;
	float mu_ = 0.7f;
	bool is_collision_;
	float collision_num_;
	glm::vec3 collision_pos_;

	std::vector<glm::vec3> r_relative_;// form vertices to the massCenter(at beginning)

public:
	RigidBody();
	RigidBody(	
		std::vector<glm::vec3> vertices,
		std::vector<Face> faces,
		float M,
		glm::mat3 I_ref,
		glm::vec3 c0,
		glm::quat rotationQuaternion = glm::quat(1.0f, 0.0f, 0.0f, 0.0f),
		glm::vec3 Vc = glm::vec3(0.0f, 0.0f, 0.0f),
		glm::vec3 omega = glm::vec3(0.0f, 0.0f, 0.0f));

	glm::mat3 rotationMatrix();
	glm::mat3 I();
	glm::mat3 I_inv();
	glm::mat3 ToCrossProduct(glm::vec3 v);

	glm::vec3 Fc();// force
	glm::vec3 Tau();// torque
	     
	void handleCollision(const Field3<float>& solid, const float& l);

	void update(float dt, const Field3<float>& solid, const float& l);

	// impulse
	void applyImpulse(const glm::vec3 f, float dt);
	void applyTorqueImpulse(const glm::vec3 tau, float dt);	

	void movePosition(const glm::vec3& dr);
	void setVelocity(const glm::vec3& vc);	

	//
	float M();
	glm::vec3 c0();
	glm::vec3 c();
	glm::vec3 Vc();
	glm::vec3 omega();
	glm::vec3 V(const glm::vec3& pos);
	std::vector<glm::vec3> vertices();
	std::vector<Face> faces();
};

#endif
