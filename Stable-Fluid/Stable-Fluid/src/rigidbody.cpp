#include "rigidbody.h"

RigidBody::RigidBody() {}

RigidBody::RigidBody(
	std::vector<glm::vec3> vertices,
	std::vector<Face> faces,
	float M,
	glm::mat3 I_ref,
	glm::vec3 c,
	glm::quat rotationQuaternion,
	glm::vec3 Vc,
	glm::vec3 omega) : vertices_(vertices), faces_(faces), M_(M), I_ref_(I_ref), I_inv_ref_(glm::inverse(I_ref)), c0_(c), c_(c), rotationQuaternion_(rotationQuaternion), Vc_(Vc), omega_(omega) {		
	L_ = I_ref * omega;
	for (int i = 0; i < vertices.size(); i++) {
		r_relative_.push_back(vertices[i] - c);
	}
}

glm::mat3 RigidBody::rotationMatrix() {
	return glm::mat3_cast(rotationQuaternion_);
}

glm::mat3 RigidBody::I() {
	return rotationMatrix() * I_ref_ * glm::transpose(rotationMatrix());
}

glm::mat3 RigidBody::I_inv() {
	return rotationMatrix() * I_inv_ref_ * glm::transpose(rotationMatrix());
}

glm::mat3 RigidBody::ToCrossProduct(glm::vec3 v) {
	return glm::mat3(
		0, v.z, -v.y,
		-v.z, 0, v.x,
		v.y, -v.x, 0
	);
}

glm::vec3 RigidBody::Fc() {
	return glm::vec3(0.0f, -M_ * G_, 0.0f) - M_ * Gamma_ * Vc_;
}

glm::vec3 RigidBody::Tau() {
	return glm::vec3(0);
}

void RigidBody::handleCollision(const Field3<float>& solid, const float& l) {
	// detecting
	is_collision_ = false;
	collision_num_ = 0;
	collision_pos_ = glm::vec3(0);
	for (int i = 0; i < r_relative_.size(); i++) {
		if (interpolate_value((c_ + rotationMatrix() * r_relative_[i]) / l, solid) <= 0) {// inside
			glm::vec3 grad;
			interpolate_gradient(grad, (c_ + rotationMatrix() * r_relative_[i]) / l, solid);
			grad = glm::normalize(grad);
			if (glm::dot(Vc_ + glm::cross(omega_, rotationMatrix() * r_relative_[i]), grad) < 0) {// v_n < 0
				is_collision_ = true;				
				collision_num_++;
				collision_pos_ += r_relative_[i];
			}
		}
	}
	if (!is_collision_) return;
	collision_pos_ /= collision_num_;
	// handling
	glm::vec3 grad;
	interpolate_gradient(grad, (c_ + rotationMatrix() * collision_pos_) / l, solid);
	grad = glm::normalize(grad);
	glm::vec3 Vi = Vc_ + glm::cross(omega_, rotationMatrix() * collision_pos_);
	glm::vec3 Vi_n = glm::dot(grad, Vi) * grad;
	glm::vec3 Vi_t = Vi - Vi_n;

	glm::mat3 K = glm::mat3(1.0f / M_) - ToCrossProduct(rotationMatrix() * Vi) * I_inv() * ToCrossProduct(rotationMatrix() * Vi);
	glm::vec3 Vi_new = -e_ * Vi_n + glm::max(1.0f - mu_ * (1.0f + e_) * glm::length(Vi_n) / glm::length(Vi_t), 0.0f) * Vi_t;
	glm::vec3 j = glm::inverse(K) * (Vi_new - Vi);

	L_ += glm::cross(rotationMatrix() * collision_pos_, j);

	// set boundary y=0.0195
	is_collision_ = false;
	collision_num_ = 0;
	collision_pos_ = glm::vec3(0);
	for (int i = 0; i < r_relative_.size(); i++) {
		if ((c_ + rotationMatrix() * r_relative_[i]).y < 0.0195f) {// inside									
			if ((Vc_ + glm::cross(omega_, rotationMatrix() * r_relative_[i])).y < 0) {// v_y < 0
				is_collision_ = true;
				collision_num_++;
				collision_pos_ += r_relative_[i];
			}
		}
	}
	if (is_collision_) {
		collision_pos_ /= collision_num_;
		c_ += (0.0195f - (c_ + rotationMatrix() * collision_pos_).y) * glm::vec3(0.0f, 1.0f, 0.0f);
	}

	Vc_ += j / M_;
	omega_ = I_inv() * L_;

	return;
}

void RigidBody::update(float dt, const Field3<float>& solid, const float& l) {
	// update Vc_ and c_
	Vc_ += Fc() / M_ * dt;
	c_ += Vc_ * dt;
	// update rotationQuaternion_ and omega_
	L_ += Tau() * dt;
	rotationQuaternion_ += glm::quat(0, omega_ * dt * 0.5f) * rotationQuaternion_;
	omega_ = I_inv() * L_;

	handleCollision(solid, l);

	return;
}

void RigidBody::applyImpulse(const glm::vec3 f, float dt) {
	Vc_ += f / M_ * dt;
	return;
}

void RigidBody::applyTorqueImpulse(const glm::vec3 tau, float dt) {
	L_ += tau * dt;
	omega_ = I_inv() * L_;
	return;
}

void RigidBody::movePosition(const glm::vec3& dr) {
	c_ += dr;
}

void RigidBody::setVelocity(const glm::vec3& vc) {
	Vc_ = vc;
}


float RigidBody::M() {
	return M_;
}

glm::vec3 RigidBody::c0() {
	return c0_;
}

glm::vec3 RigidBody::c() {
	return c_;
}

glm::vec3 RigidBody::Vc() {
	return Vc_;
}

glm::vec3 RigidBody::omega() {
	return omega_;
}

glm::vec3 RigidBody::V(const glm::vec3& pos) {
	return Vc_ + glm::cross(omega_, pos - c_);
}

std::vector<glm::vec3> RigidBody::vertices() {
	std::vector<glm::vec3> vertices;
	for (int i = 0; i < r_relative_.size(); i++) {
		vertices.push_back(c_ + rotationMatrix() * r_relative_[i]);
	}
	return vertices;
}

std::vector<Face> RigidBody::faces() {
	return faces_;
}
