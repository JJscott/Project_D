
#include <GLFW/glfw3.h>

#include <cstring>
#include <list>
#include <map>
#include <algorithm>

#include "glfw_helper.h"

using namespace std;

namespace glfwpp {
	
	class keyboard_state {
	private:
		bool m_down[GLFW_KEY_LAST];
		
	public:
		void clear() {
			memset(m_down, 0, sizeof(bool) * GLFW_KEY_LAST);
		}
		
		keyboard_state() {
			clear();
		}
		
		void press(int key) {
			m_down[key] = true;
		}
		
		void release(int key) {
			m_down[key] = false;
		}
		
		bool down(int key) {
			return m_down[key];
		}
	};
	
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
	
	map<GLFWwindow *, keyboard_state *> kb_states;
	
	list<key_listener_binding> kl_bindings;
	list<mouse_listener_binding> ml_bindings;
	
	void handle_focus(GLFWwindow *window, int focused) {
		// clear key state if unfocused
		if (!focused) {
			kb_states[window]->clear();
		}
	}
	
	void handle_key(GLFWwindow *window, int key, int scancode, int action, int mods) {
		for (list<key_listener_binding>::const_iterator it = kl_bindings.begin(); it != kl_bindings.end(); it++) {
			if (it->window == NULL || it->window == window) {
				it->listener->handleKey(window, key, scancode, action, mods);
			}
		}
		// update keyboard state
		if (key != GLFW_KEY_UNKNOWN) {
			if (action == GLFW_PRESS || action == GLFW_REPEAT) kb_states[window]->press(key);
			if (action == GLFW_RELEASE) kb_states[window]->release(key);
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
	void setCallbacks(GLFWwindow *window) {
		// input callbacks
		glfwSetKeyCallback(window, handle_key);
		glfwSetCharCallback(window, handle_char);
		glfwSetMouseButtonCallback(window, handle_mouse_button);
		glfwSetCursorPosCallback(window, handle_cursor_pos);
		glfwSetCursorEnterCallback(window, handle_cursor_enter);
		glfwSetScrollCallback(window, handle_scroll);
		// focus callback
		glfwSetWindowFocusCallback(window, handle_focus);
		// setup the keyboard state
		if (kb_states.find(window) == kb_states.end()) {
			kb_states[window] = new keyboard_state();
		}
	}
	
	// check the status of a key
	bool getKey(GLFWwindow * window, int key) {
		if (kb_states.find(window) == kb_states.end()) return false;
		return kb_states[window]->down(key);
	}
	
	// check the status of a key, then 'release' it
	bool pollKey(GLFWwindow *window, int key) {
		if (kb_states.find(window) == kb_states.end()) return false;
		bool down = kb_states[window]->down(key);
		kb_states[window]->release(key);
		return down;
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
		for (list<key_listener_binding>::iterator it = kl_bindings.begin(); it != kl_bindings.end(); it++) {
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
		for (list<mouse_listener_binding>::iterator it = ml_bindings.begin(); it != ml_bindings.end(); it++) {
			if (*it == mlb || (window == NULL && it->listener == ml)) {
				ml_bindings.erase(it);
				break;
			}
		}
	}
}
