/*
 * Camera Control
 *
 *
 */

#ifndef CAMERA_H
#define CAMERA_H

#include <GLFW/glfw3.h>

#include "glfw_helper.h"
#include "initial3d.h"

class Camera {
public:
	virtual void update() = 0;
	virtual initial3d::mat4d getViewMatrix() = 0;
	virtual ~Camera() { }
};

class FPSCamera : private initial3d::Uncopyable, public Camera, protected glfwpp::KeyListener {
private:
	GLFWwindow *window;
	initial3d::vec3d pos;
	double rot_h, rot_v;
	bool mouse_captured;
	
public:
	FPSCamera(GLFWwindow *window_, const initial3d::vec3d &pos_, double rot_h_ = 0, double rot_v_ = 0) :
		window(window_), pos(pos_), rot_h(rot_h_), rot_v(rot_v_), mouse_captured(false) {
			glfwpp::addKeyListener(window, this);
	}
	
	initial3d::quatd getOrientation() {
		return initial3d::quatd::axisangle(initial3d::vec3d::j(), rot_h) * initial3d::quatd::axisangle(initial3d::vec3d::i(), rot_v);
	}
	
	virtual void update() {
		// pixels per 2*pi
		double rot_speed = 600;
		
		if (mouse_captured) {
			int w, h;
			glfwGetWindowSize(window, &w, &h);
			double x, y;
			glfwGetCursorPos(window, &x, &y);
			x -= w * 0.5;
			y -= h * 0.5;
			rot_h += -x / rot_speed;
			rot_v += -y / rot_speed;
			rot_v = initial3d::math::clamp(rot_v, -0.499 * initial3d::math::pi(), 0.499 * initial3d::math::pi());
			glfwSetCursorPos(window, w * 0.5, h * 0.5);
		}
		
	}
	
	virtual initial3d::mat4d getViewMatrix() {
		return initial3d::mat4d::rotate(!getOrientation()) * initial3d::mat4d::translate(-pos);
	}
	
	virtual void handleKey(GLFWwindow *window, int key, int scancode, int action, int mods) {
		if (key == GLFW_KEY_GRAVE_ACCENT && action == GLFW_PRESS) mouse_captured = !mouse_captured;
	}

	virtual ~FPSCamera() {
		glfwpp::removeKeyListener(NULL, this);
	}
};

#endif // CAMERA_H
