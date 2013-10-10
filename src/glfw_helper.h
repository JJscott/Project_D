
#ifndef GLFW_HELPER_H
#define GLFW_HELPER_H

#include <GLFW/glfw3.h>

namespace glfwpp {
	
	class KeyListener {
	public:
		virtual void handleKey(GLFWwindow *window, int key, int scancode, int action, int mods) { }
		virtual void handleChar(GLFWwindow *window, unsigned int character) { }
		virtual ~KeyListener() { }
	};
	
	class MouseListener {
	public:
		virtual void handleMouseButton(GLFWwindow *window, int button, int action, int mods) { }
		virtual void handleCursorPos(GLFWwindow *window, double xpos, double ypos) { }
		virtual void handleCursorEnter(GLFWwindow *window, int entered) { }
		virtual void handleScroll(GLFWwindow *window, double xoffset, double yoffset) { }
		virtual ~MouseListener() { }
	};
	
	// set the GLFW callbacks to the ones required to use these helpers
	// uses key + mouse + focus
	// will not be called automatically
	void setCallbacks(GLFWwindow *window);
	
	// if window == NULL, listens to all windows (that the callbacks are set for)
	void addKeyListener(GLFWwindow *window, KeyListener *kl);
	
	// if window == NULL, remove from all windows
	void removeKeyListener(GLFWwindow *window, KeyListener *kl);
	
	// if window == NULL, listens to all windows (that the callbacks are set for)
	void addMouseListener(GLFWwindow *window, MouseListener *ml);
	
	// if window == NULL, remove from all windows
	void removeMouseListener(GLFWwindow *window, MouseListener *ml);
}

#endif
