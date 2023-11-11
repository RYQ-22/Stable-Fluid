#include "camera.h"

Camera::Camera(glm::vec3 position, glm::vec3 up, float yaw, float pitch) {
	front_ = glm::vec3(0.0f, 0.0f, -1.0f);
	movementSpeed_ = SPEED;
	mouseSensitivity_ = SENSITIVITY;
	zoom_ = ZOOM;
	position_ = position;
	worldUp_ = up;
	yaw_ = yaw;
	pitch_ = pitch;
	updateCameraVectors();
}

Camera::Camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw, float pitch) {
	front_ = glm::vec3(0.0f, 0.0f, -1.0f);
	movementSpeed_ = SPEED;
	mouseSensitivity_ = SENSITIVITY;
	zoom_ = ZOOM;
	position_ = glm::vec3(posX, posY, posZ);
	worldUp_ = glm::vec3(upX, upY, upZ);
	yaw_ = yaw;
	pitch_ = pitch;
	updateCameraVectors();
}

glm::mat4 Camera::GetViewMatrix() {
	return glm::lookAt(position_, position_ + front_, up_);
}

void Camera::ProcessKeyboard(Camera_Movement direction, float deltaTime) {
	float velocity = movementSpeed_ * deltaTime;
	if (direction == FORWARD)
		position_ += front_ * velocity;
	if (direction == BACKWARD)
		position_ -= front_ * velocity;
	if (direction == LEFT)
		position_ -= right_ * velocity;
	if (direction == RIGHT)
		position_ += right_ * velocity;
}

void Camera::ProcessMouseMovement(float xoffset, float yoffset, GLboolean constrainPitch) {
	xoffset *= mouseSensitivity_;
	yoffset *= mouseSensitivity_;
	yaw_ += xoffset;
	pitch_ -= yoffset;

	if (constrainPitch) {
		if (pitch_ > 89.0f)
			pitch_ = 89.0f;
		if (pitch_ < -89.0f)
			pitch_ = -89.0f;
	}
	updateCameraVectors();
}

void Camera::updateCameraVectors() {
	glm::vec3 front;
	front.x = cos(glm::radians(yaw_)) * cos(glm::radians(pitch_));
	front.y = sin(glm::radians(pitch_));
	front.z = sin(glm::radians(yaw_)) * cos(glm::radians(pitch_));
	front_ = glm::normalize(front);
	right_ = glm::normalize(glm::cross(front_, worldUp_));
	up_ = glm::normalize(glm::cross(right_, front_));
}
