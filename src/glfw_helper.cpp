
#include <GLFW/glfw3.h>

#include <list>
#include <algorithm>

#include "glfw_helper.h"

using namespace std;

namespace glfwpp {
	
	struct key_listener_binding {
		GLFWwindow *window;
		KeyListener *listener;
		
		bool operator==(const key_listener_binding &rhs) const {
			return window == rhs.window && listener == rhs.listener;
		}
	};

	struct mouse_listener_binding {
		GLFWwindow *window;
		MouseListener *listener;

		bool operator==(const mouse_listener_binding &rhs) const {
			return window == rhs.window && listener == rhs.listener;
		}
	};

	list<key_listener_binding> kl_bindings;
	list<mouse_listener_binding> ml_bindings;

	void handle_key(GLFWwindow *window, int key, int scancode, int action, int mods) {
		for (list<key_listener_binding>::const_iterator it = kl_bindings.begin(); it != kl_bindings.end(); it++) {
			if (it->window == NULL || it->window == window) {
				it->listener->handleKey(window, key, scancode, action, mods);
			}
		}
	}

	void handle_char(GLFWwindow *window, unsigned int character) {
		for (list<key_listener_binding>::const_iterator it = kl_bindings.begin(); it != kl_bindings.end(); it++) {
			if (it->window == NULL || it->window == window) {
				it->listener->handleChar(window, character);
			}
		}
	}

	void handle_mouse_button(GLFWwindow *window, int button, int action, int mods) {
		for (list<mouse_listener_binding>::const_iterator it = ml_bindings.begin(); it != ml_bindings.end(); it++) {
			if (it->window == NULL || it->window == window) {
				it->listener->handleMouseButton(window, button, action, mods);
			}
		}
	}

	void handle_cursor_pos(GLFWwindow *window, double xpos, double ypos) {
		for (list<mouse_listener_binding>::const_iterator it = ml_bindings.begin(); it != ml_bindings.end(); it++) {
			if (it->window == NULL || it->window == window) {
				it->listener->handleCursorPos(window, xpos, ypos);
			}
		}
	}

	void handle_cursor_enter(GLFWwindow *window, int entered) {
		for (list<mouse_listener_binding>::const_iterator it = ml_bindings.begin(); it != ml_bindings.end(); it++) {
			if (it->window == NULL || it->window == window) {
				it->listener->handleCursorEnter(window, entered);
			}
		}
	}

	void handle_scroll(GLFWwindow *window, double xoffset, double yoffset) {
		for (list<mouse_listener_binding>::const_iterator it = ml_bindings.begin(); it != ml_bindings.end(); it++) {
			if (it->window == NULL || it->window == window) {
				it->listener->handleScroll(window, xoffset, yoffset);
			}
		}
	}

	// set the GLFW callbacks to the ones required to use these helpers
	// will not be called automatically
	void setInputCallbacks(GLFWwindow *window) {
		glfwSetKeyCallback(window, handle_key);
		glfwSetCharCallback(window, handle_char);
		glfwSetMouseButtonCallback(window, handle_mouse_button);
		glfwSetCursorPosCallback(window, handle_cursor_pos);
		glfwSetCursorEnterCallback(window, handle_cursor_enter);
		glfwSetScrollCallback(window, handle_scroll);
	}
	
	// if window == NULL, listens to all windows
	void addKeyListener(GLFWwindow *window, KeyListener *kl) {
		key_listener_binding klb;
		klb.window = window;
		klb.listener = kl;
		if (find(kl_bindings.begin(), kl_bindings.end(), klb) == kl_bindings.end()) {
			// dont add duplicates
			kl_bindings.push_back(klb);
		}
	}
	
	// if window == NULL, remove from all windows
	void removeKeyListener(GLFWwindow *window, KeyListener *kl) {
		key_listener_binding klb;
		klb.window = window;
		klb.listener = kl;
		for (list<key_listener_binding>::const_iterator it = kl_bindings.begin(); it != kl_bindings.end(); it++) {
			if (*it == klb || (window == NULL && it->listener == kl)) {
				kl_bindings.erase(it);
				break;
			}
		}
	}
	
	// if window == NULL, listens to all windows
	void addMouseListener(GLFWwindow *window, MouseListener *ml) {
		mouse_listener_binding mlb;
		mlb.window = window;
		mlb.listener = ml;
		if (find(ml_bindings.begin(), ml_bindings.end(), mlb) == ml_bindings.end()) {
			// dont add duplicates
			ml_bindings.push_back(mlb);
		}
	}
	
	// if window == NULL, remove from all windows
	void removeMouseListener(GLFWwindow *window, MouseListener *ml) {
		mouse_listener_binding mlb;
		mlb.window = window;
		mlb.listener = ml;
		for (list<mouse_listener_binding>::const_iterator it = ml_bindings.begin(); it != ml_bindings.end(); it++) {
			if (*it == mlb || (window == NULL && it->listener == ml)) {
				ml_bindings.erase(it);
				break;
			}
		}
	}
}
