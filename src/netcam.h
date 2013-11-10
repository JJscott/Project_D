
#include <iostream>

#include "initial3d.h"
#include "atmos.h"
#include "camera.h"
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef NETCAM_H
#define NETCAM_H


class NetworkCamera : public Camera {
private:
	initial3d::mat4d m_cammatrix;


public:
	NetworkCamera() {
		using namespace initial3d;
		m_cammatrix = mat4d::translate(vec3d(0, atmos::Rg + 1000, 0));
		std::cout << m_cammatrix << std::endl;
	}

	virtual void update() {
		
	}

	virtual initial3d::mat4d getViewMatrix() {
	  static int setup=0;
	  static int sock=0;
	  static socklen_t length=0;
	  static struct sockaddr_in server;
	  static struct sockaddr_in client;
	  char buf[1024];

	  if(setup==0){
	    /* Create socket from which to read. */

	    sock = socket(AF_INET, SOCK_DGRAM, 0);
	    if (sock < 0) {
	      perror("opening datagram socket");
	      exit(1);
	    }
	    /* Create server with wildcards. */
	    server.sin_family = AF_INET;
	    server.sin_addr.s_addr = INADDR_ANY;
	    server.sin_port = 0;

	    if (bind(sock, (struct sockaddr*)&server, sizeof(server))) {
	      perror("binding datagram socket");
	      exit(1);
	    }

	    /* Find assigned port value and print it out. */
	    length = sizeof(server);
	    if (getsockname(sock, (struct sockaddr*)&server, &length)) {
	      perror("getting socket server");
	      exit(1);
	    }

	    fprintf(stderr,"Socket has port #%d\n", ntohs(server.sin_port));
	    setup=1;
	  }
                 
		 

	  /* Read from the socket */
	  length = sizeof(client);
	  if (recvfrom(sock, buf, 1024, 0, (struct sockaddr*)&client, &length) < 0) perror("receiving datagram packet");
	  else 
	    {
	      //TODO divide floats by 63.0 and place in matrix 63.0 + me = 64.0
	      printf("Client says: -->%f %f %f %f %f %f %f %f %f\n", ((float)buf[0])/63.0, ((float)buf[4])/63.0, ((float)buf[8])/63.0, ((float)buf[12])/63.0, ((float)buf[16])/63.0, ((float)buf[20])/63.0, ((float)buf[24])/63.0, ((float)buf[28])/63.0, ((float)buf[32])/63.0);
	    }
	  /* Send message. */
	  //sprintf(buf, "%d", (atoi(buf))+1);
	  //if (sendto(sock, buf, sizeof(buf), 0, (struct sockaddr*)&client, sizeof(client)) < 0)
	  //perror("sending datagram message");

		
	  //close(sock);
       
		return !m_cammatrix;
	}

};


#endif
