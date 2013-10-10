/*
 * Camera Control
 *
 *
 */

#ifndef CAMERA_H
#define CAMERA_H

#include <GLFW/glfw3.h>

#include <iostream>

#include "glfw_helper.h"
#include "initial3d.h"
#include "hrtime.h"

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
	initial3d::quatd ori;
	double rot_v;
	bool mouse_captured;
	hrtime_t time_last;
	double speed;
	
public:
	FPSCamera(GLFWwindow *window_, const initial3d::vec3d &pos_, double rot_h_ = 0, double rot_v_ = 0) :
		window(window_), pos(pos_), ori(), rot_v(rot_v_), mouse_captured(false) {
			glfwpp::addKeyListener(window, this);
			getHighResTime(&time_last);
			speed = 1000;
	}
	
	initial3d::quatd getOrientation() {
		return ori * initial3d::quatd::axisangle(initial3d::vec3d::i(), rot_v);
	}
	
	virtual void update() {
		using namespace initial3d;
		using namespace std;
		
		// time since last update
		hrtime_t time_now;
		getHighResTime(&time_now);
		double dt = highResTimeToSec(&time_last, &time_now);
		time_last = time_now;

		// pixels per 2*pi
		double rot_speed = 600;
		
		vec3d up = ~pos; // we won't be going to the centre of the planet, right?
		vec3d forward = -~(ori * vec3d::k()).reject(up);
		vec3d side = ~(forward ^ up);

		if (mouse_captured) {
			int w, h;
			glfwGetWindowSize(window, &w, &h);
			double x, y;
			glfwGetCursorPos(window, &x, &y);
			x -= w * 0.5;
			y -= h * 0.5;
			double rot_h = -x / rot_speed;
			ori = quatd::axisangle(up, rot_h) * ori;
			rot_v += -y / rot_speed;
			rot_v = math::clamp(rot_v, -0.499 * math::pi(), 0.499 * math::pi());
			glfwSetCursorPos(window, w * 0.5, h * 0.5);
		}
		
		vec3d move = vec3d::zero();

		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) move += forward;
		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) move -= forward;
		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) move -= side;
		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) move += side;
		if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) move -= up;
		if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) move += up;

		try {
			double pos_mag = +pos;
			vec3d dpos = ~move * speed * dt;
			vec3d pos2 = pos + dpos;
			try {
				quatd rot = quatd::axisangle(~pos ^ ~pos2, pos.angle(pos2));
				ori = rot * ori;
			} catch (nan_error &e) {
				// no rotation, do nothing
			}
			pos = pos2;
			+pos = pos_mag + dpos * up;
			
			cout << pos << endl;
			
		} catch (nan_error &e) {
			// no movement, do nothing
		}

	}
	
	virtual initial3d::mat4d getViewMatrix() {
		return initial3d::mat4d::rotate(!getOrientation()) * initial3d::mat4d::translate(-pos);
	}
	
	virtual void handleKey(GLFWwindow *window, int key, int scancode, int action, int mods) {
		if (action == GLFW_PRESS) {
			if (key == GLFW_KEY_GRAVE_ACCENT) mouse_captured = !mouse_captured;
			if (key == GLFW_KEY_LEFT_BRACKET) speed *= 0.5;
			if (key == GLFW_KEY_RIGHT_BRACKET) speed *= 2;
		}
	}

	virtual ~FPSCamera() {
		glfwpp::removeKeyListener(NULL, this);
	}
};

#endif // CAMERA_H
