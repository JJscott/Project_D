
#include <GLFW/glfw3.h>

#include "glfw_helper.h"

namespace glfwpp {
	
	// set the GLFW callbacks to the ones required to use these helpers
	// will not be called automatically
	void init() {
		
	}
	
	// if window == NULL, listens to all windows
	void addKeyListener(GLFWwindow *window, KeyListener *kl) {
		
	}
	
	// if window == NULL, remove from all windows
	void removeKeyListener(GLFWwindow *window, KeyListener *kl) {
		
	}
	
	// if window == NULL, listens to all windows
	void addMouseListener(GLFWwindow *window, MouseListener *ml) {
		
	}
	
	// if window == NULL, remove from all windows
	void removeMouseListener(GLFWwindow *window, MouseListener *ml) {
		
	}
}
