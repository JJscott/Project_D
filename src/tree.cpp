
#include <vector>
#include <set>
#include <iostream>
#include <map>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include "GLee.h"
#include <GL/glu.h>
#include "tree.h"
#include "volume.h"
#include "initial3d.h"

using namespace std;
using namespace vol;
using namespace initial3d;


// k subgraph resolution (number of nodes)
// d node linking distance
// a subgraph shrinking parameter (a < 1)
// A subgraph sizing angle
// Note, a node N has property~x, spatial position.



// 	V (graph volume),
// 	R (root node),
// 	b (branching factor),
// 	K (number of nodes in graph),
// 	r (size),
// 	L (remaining lifespan)










TreeGenerator::TreeGenerator() {
	//set global variables here

	param_nodeLinkDistance = 0.15;
	param_subgraphSize = 0.5;
	param_subgraphAngle = 0.0;

	m_maxDepth = 3;
	m_branchFactor = 5;
	m_graphPopulation = 10000;

	m_rand = new initial3d::random(33);
}


// TreeGenerator(graphResolution, linkingDistance, shrinkingParameter, sizingAngle) {
// 	//set global variables here
// }


// //creates a mesh for the given tree and maps it to a display list given by the GLuint id
// void TreeGenerator::buildTree(GLuint id) {


// 	// TreeNode *root = new TreeNode(vec3d(0,0,0));
// 	// TreeNode *child = new TreeNode(vec3d(0,4,0));
// 	// TreeNode *child1 = new TreeNode(vec3d(0,6,-1));
// 	// TreeNode *child2 = new TreeNode(vec3d(0.5,6,1));
// 	// TreeNode *child3 = new TreeNode(vec3d(0.5,6,-0.5));

// 	// root->children.push_back(child);
// 	// child->children.push_back(child1);
// 	// child->children.push_back(child2);
// 	// child->children.push_back(child3);

// 	// root->calculateBranchCount();
// 	// root->calculateRadius(0.1, 2.0);
// 	// createMesh(id, root);


// 	// Volume *s1 = new Sphere(vec3d::j() * 2.0, 2.0);
// 	// Volume *s2 = new Sphere(vec3d::j() * -1.0, 1.0);
// 	// Volume *volume = new Union(s1, s2);



// 	// Volume *volume = new Sphere(vec3d::zero(), 2.0);
// 	// Volume *volume = new Cone(vec3d(0,0,-5.0), 5.0, 5.0);
// 	// Volume *volume = new Cone(vec3d::zero(), 2.0, 1.0);
// 	// Volume *volume = new Cylinder(vec3d(0,0,-2.0), 2.0, 1.0);
// 	// Volume *volume = new Cylinder(vec3d::zero(), 2.0, 1.0);
// 	// Volume *volume = new Box(vec3d(-2,-2,-2), vec3d(2,2,2));
// 	// Volume *volume = new Sphere(vec3d::j(), 2.0);
// 	// Volume *volume = new Sphere(vec3d::j(), 2.0);
// 	// Volume *volume = new Sphere(vec3d::j(), 2.0);
// 	// Volume *volume = new Sphere(vec3d::j(), 2.0);


// 	quatd rot(vec3d::i(), -math::pi()/2);
// 	Volume *cyl = new Cylinder(vec3d::zero(), 6.0, 0.25);
// 	// Volume *volume = new Rotation(cyl, rot);
// 	Volume *rotate = new Rotation(cyl, rot);
// 	Volume *sphere = new Sphere(vec3d(0, 5, 0), 2.0);
// 	// Volume *sphere = new Sphere(vec3d(0, 2, 0), 2.0);
// 	Volume *volume = new Union(rotate, sphere);

// 	// Volume *volume = new Rotation(uni, rot, vec3d(0, 6, 0));

// 	cout << "VMin " << volume->getMin() << " : VMax " << volume->getMax() << endl;



// 	Node *root = new Node(vec3d(0.0,0.0,0.0));
// 	Graph *graph = new Graph();
// 	graph->addAndLinkNode(root, 0);
// 	map<Node *, Node *> parents;

// 	// createGraph(volume, graph);
// 	createGraph2(volume, graph);
// 	// drawGraphPoints(id, graph);


// //--------------------------------------------------------------------

// 	mapParents(graph, root, parents);

// 	// create some end nodes
// 	vector<Node *> endPoints;
// 	int factor = min(m_branchFactor, (int)graph->nodes.size());
// 	while (factor > 0) {
// 		int i = m_rand->nextInt(graph->nodes.size());
// 		if (graph->nodes[i] != root && find(endPoints.begin(), endPoints.end(), graph->nodes[i]) == endPoints.end()) {
// 			endPoints.push_back(graph->nodes[i]);
// 			factor--;
// 		}
// 	}

// 	// TreeNode *t_root = buildNaiveTree(graph, root, parents);
// 	TreeNode *t_root = buildTree(graph, root, parents, endPoints);


// 	t_root->calculateBranchCount();
// 	t_root->calculateRadius(0.1, 2.0);
// 	createMesh(id, t_root);
// 	// cout << graph->nodes.size() <<endl;

// //--------------------------------------------------------------------

// }


void TreeGenerator::buildTree(GLuint id) {
	quatd rot = quatd::axisangle(vec3d::i(), -math::pi()/2);
	Volume *cyl = new Cylinder(vec3d::zero(), 3, 0.25);
	Volume *rotate = new Rotation(cyl, rot);
	Volume *sphere = new Sphere(vec3d(0, 3, 0), 1.0);
	Volume *volume = new Union(rotate, sphere);


	Node *root = new Node(vec3d(0.0,0.0,0.0));
	Graph *graph = new Graph();
	graph->addAndLinkNode(root, 0);
	map<Node *, Node *> parents;

	cout << "Creating graph" << endl;

	createGraph2(volume, graph);
	mapParents(graph, root, parents);

	cout << "Soring out shit" << endl;

	// create some end nodes
	vector<Node *> endPoints;
	int factor = min(m_branchFactor, (int)graph->nodes.size());
	while (factor > 0) {
		int i = m_rand->nextInt(graph->nodes.size());
		if (graph->nodes[i] != root && find(endPoints.begin(), endPoints.end(), graph->nodes[i]) == endPoints.end()) {
			endPoints.push_back(graph->nodes[i]);
			factor--;
		}
	}

		// TreeNode *t_root = buildNaiveTree(graph, root, parents);
	TempTreeNode *temp_root = new TempTreeNode(root);
	vector<TempTreeNode *> endPointContainers;

	buildTree(graph, temp_root, parents, endPoints, endPointContainers);

	cout << "Begin the PAIN" << endl;
	buildSubTree(temp_root, endPointContainers, 1);
	

	TreeNode *t_root = temp_root->toTreeNode();


	t_root->calculateBranchCount();
	t_root->calculateRadius(0.15, 2.0);
	createMesh(id, t_root);
}





void TreeGenerator::buildSubTree(TempTreeNode *parent_root, vector<TempTreeNode *> &endPointContainers, int depth) {

	cout << "depth : " << depth << " : size : " << endPointContainers.size() << endl;

	for (int i = 0; i < (int) endPointContainers.size(); i++) {
		cout << "    :: count " << i << endl;

		Graph *subGraph = new Graph();
		Node *subRoot = endPointContainers[i]->node;

		vec3d ori = ~(subRoot->position - parent_root->node->position);
		
		double dot = -vec3d::k() * ori;
		vec3d crossp = -vec3d::k() ^ ori;
		double x = crossp.x();
		double y = crossp.y();
		double z = crossp.z();

		// quatd rot1(dot, x, y, z);
		quatd rot1;
		try {
			rot1 = quatd::axisangle(crossp, math::acos(dot));

		} catch (nan_error &e) {}

		double scale = math::pow<double>(param_subgraphSize, depth);

		Volume *cone = new Cone(subRoot->position - vec3d::k(2 * scale), 2 *scale, 3 * scale);
		Volume *volume = new Rotation(cone, rot1, subRoot->position);

		subGraph->addAndLinkNode(subRoot, 0);
		map<Node *, Node *> parents;


		createGraph2(volume, subGraph);
		mapParents(subGraph, subRoot, parents);

		// create some end nodes
		vector<Node *> endPoints;
		vector<TempTreeNode *> newEndPointContainers;
		int factor = min(m_branchFactor, (int)subGraph->nodes.size()-1);
		while (factor > 0) {
			int j = m_rand->nextInt(subGraph->nodes.size());
			if (subGraph->nodes[j] != subRoot && find(endPoints.begin(), endPoints.end(), subGraph->nodes[j]) == endPoints.end()) {
				endPoints.push_back(subGraph->nodes[j]);
				factor--;
			}
		}
		buildTree(subGraph, endPointContainers[i], parents, endPoints, newEndPointContainers);

		if (depth < m_maxDepth) {
			for (int j = 0; j<(int)newEndPointContainers.size(); j++) {
				buildSubTree(endPointContainers[i], newEndPointContainers, depth+1);
			}
		}
		cout << "    :: again " << i << endl;
	}

}















void TreeGenerator::createGraph(Volume *volume, Graph *graph) {

	vec3d min = volume->getMin();
	vec3d max = volume->getMax();
	vec3d scale = max - min;
	if (math::isinf(min.x() + min.y() + min.z()) || math::isinf(max.x() + max.y() + max.z()))
		return; //or error out

	// int estNumPoints = (scale.x() * scale.y() * scale.z());
	// estNumPoints = 8 * estNumPoints / param_nodeLinkDistance;
 
	int estNumPoints = 100000;

	for (int i = 0; i < estNumPoints; i++) {

		double x = m_rand->nextDouble() * scale.x() + min.x();
		double y = m_rand->nextDouble() * scale.y() + min.y();
		double z = m_rand->nextDouble() * scale.z() + min.z();
		vec3d point(x,y,z);
		bool reject = false;

		if (volume->contains(point)) {
			for (int i=0; i<(int)graph->nodes.size(); i++) {
				if (+(point - graph->nodes[i]->position) < param_nodeLinkDistance) {
					reject = true;
					continue;
				}
			}

			if (!reject) {
				graph->addAndLinkNode(new Node(point), param_nodeLinkDistance * 2);
			}
		}
	}

	// for (int i=0; i<(int)graph->nodes.size(); i++) {
	// 	cout << "Node id " << graph->nodes[i] << " : pos " << graph->nodes[i]->position << endl;
	// }

	// for (int i=0; i<(int)graph->edges.size(); i++) {
	// 	cout << "Edge id " << graph->edges[i] << " : n1 " << graph->edges[i]->n1 << " : n2 " << graph->edges[i]->n2 << endl;
	// }
}




void TreeGenerator::createGraph2(Volume *volume, Graph *graph) {
	VolumeGrid(volume, param_nodeLinkDistance, m_rand).populate(graph);
}


//
// I DONT KNOW HOW THE F*** THIS WORKS BUT IT SORT OF DOES
//
// HOWEVER IT DOES NOT RUN IN VS DEBUG MODE
//
void TreeGenerator::mapParents(Graph *graph, Node *root, map<Node *, Node *> &parent) {
	priority_queue<NodeP> allNodes;
	map<Node*, double> cost;

	//create NodeP for all nodes in the graph
	for (int i = 0; i<(int)graph->nodes.size(); i++){
		if (graph->nodes[i] != root) {
			cost[graph->nodes[i]] = math::inf<double>();
		}
		else {
			cost[graph->nodes[i]] = 0.0;
		}
		allNodes.push(NodeP(graph->nodes[i], &cost));
	}

	while (allNodes.size() > 0) {
		NodeP next = allNodes.top();
		allNodes.pop();

		if (math::isinf(cost[next.node]))
			continue;

		for (int i = 0; i<(int)next.node->edges.size(); i++) {

			Edge * edge = next.node->edges[i];
			double alt = cost[next.node] + edge->weight;
			Node *neighbour = edge->getOther(next.node);

			if (alt < cost[neighbour]) {
				cost[neighbour] = alt;
				parent[neighbour] = next.node;
			}
		}
	}
}

//
// THIS ONE WORKS PROPERLY BUT IS REALLY SLOW
//
//void TreeGenerator::mapParents(Graph *graph, Node *root, map<Node *, Node *> &parent) {
//	priority_queue<NodeP> fringe;
//	unordered_map<Node *, double> cost;
//
//	// set root to cost 0, add to fringe
//	cost[root] = 0.0;
//	fringe.push(NodeP(root, 0));
//
//	while (fringe.size() > 0) {
//		NodeP next = fringe.top();
//		fringe.pop();
//
//		// lower cost already found for this node, ignore
//		if (cost[next.node] < next.cost) continue;
//
//		for (int i=0; i<(int)next.node->edges.size(); i++) {
//			Edge *edge = next.node->edges[i];
//			Node *neighbour = edge->getOther(next.node);
//			auto cost_it = cost.find(neighbour);
//			// current best cost to neighbour, inf if not visited yet
//			double cost0 = cost_it == cost.end() ? math::inf<double>() : cost_it->second;
//			// cost to neighbour through this node
//			double cost1 = cost[next.node] + edge->weight;
//
//			if (cost1 < cost0) {
//				// cost through this node better than current best
//				cost[neighbour] = cost1;
//				parent[neighbour] = next.node;
//				// add neighbour to fringe (some fringe elements may still refer to this node with a higher cost, they will be ignored)
//				fringe.push(NodeP(neighbour, cost1));
//			}
//		}
//	}
//}


// Tree * mapEndPoints(Graph *graph, Node *root, map<Node *, Node *> &parent) {
// 	// get rid of all nodes without parents (exlcuding root)
// 	// for (int i=0; i<(int)graph->nodes().size(); i++) {
// 	// 	Node *node = graph->nodes()[i];
// 	// 	if (node != root && parent.find(node) == parent.end()) {
// 	// 		graph.erase();
// 	// 	}
// 	// }
// 	// m_rand->nextDouble();
// }


void TreeGenerator::buildTree(Graph *graph, TempTreeNode *tree, map<Node *, Node *> &parent, vector<Node*> &endPoints, vector<TempTreeNode *> &endPointContainers) {
	map<Node *, TempTreeNode *> container;
	container[tree->node] = tree;

	for (int i=0; i<(int)endPoints.size(); i++) {
		Node *child_n = endPoints[i];
		TempTreeNode *child_t;
		Node * parent_n = parent[child_n];
		TempTreeNode *parent_t;

		//there exists no container for this end node
		if (container.find(child_n) == container.end()) {
			child_t = new TempTreeNode(child_n);
			endPointContainers.push_back(child_t);
			container[child_n] = child_t;

			//until we reach the existing container path
			//create new containers along the way
			while (container.find(parent_n) == container.end() && parent_n !=tree->node) {
				parent_t = new TempTreeNode(parent_n);
				container[parent_n] = parent_t;
				parent_t->children.push_back(child_t);
				child_n = parent_n;
				child_t = parent_t;
				parent_n = parent[child_n];
			}

			//reached the container line
			//and the current attachement to the line
			parent_t = container[parent_n];
			parent_t->children.push_back(child_t);
		} else {
			endPointContainers.push_back(container[child_n]);
		}
	}
}


TreeNode * TreeGenerator::buildNaiveTree(Graph *graph, Node *root, map<Node *, Node *> &parent) {
	map<Node *, TempTreeNode *> container;
	TempTreeNode *tree = new TempTreeNode(root);
	container[root] = tree;

	for (int i=0; i<(int)graph->nodes.size(); i++) {
		Node *node = graph->nodes[i];
		//for all connected nodes ^_^
		if (node != root && parent.find(node) != parent.end()) {
			Node * nParent = parent[node];
			TempTreeNode *tNode;
			//if already exists
			if (container.find(node) != container.end())
				tNode = container[node];
			else {
				tNode = new TempTreeNode(node);
				container[node] = tNode;
			}


			if (container.find(nParent) != container.end()) //parent exsists, add
				container[nParent]->children.push_back(tNode);
			else {//else create and stuff
				TempTreeNode *tParent = new TempTreeNode(nParent);
				container[nParent] = tParent;
				tParent->children.push_back(tNode);
			}
		}
	}
	return tree->toTreeNode();
}


void TreeGenerator::createMesh(GLuint id, TreeNode *root) {
	glNewList(id, GL_COMPILE);
	GLUquadric *quad = gluNewQuadric();
	glPushMatrix();
		drawBranch(quad, root);
	glPopMatrix();
	gluDeleteQuadric(quad);
	glEndList();
}

void TreeGenerator::drawBranch(GLUquadric *quad, TreeNode *root) {
	// cout << " pos " << root->position << " : radius " << root->radius << " : count " << root->branchCount << endl;
	for (int i = 0; i < (int)root->children.size(); i++) {
		drawCylinder(quad, root->position, root->children[i]->radius, root->children[i]->position, root->children[i]->radius);
		drawBranch(quad, root->children[i]);
	}
}

void TreeGenerator::drawCylinder(GLUquadric *quad, const vec3d &from, double base, const vec3d &to, double top) {
	vec3d final = ~(to - from);
	double w = acos(vec3d::k() * final) * 180.0 / math::pi();
	vec3d cp = vec3d::k() ^ final;

	glPushMatrix();
		glTranslated(from.x(), from.y(), from.z());
		glRotated(w, cp.x(), cp.y(), cp.z());
		gluCylinder(quad, base, top, +(from - to), 20, 20);
	glPopMatrix();
}

void TreeGenerator::drawGraphPoints(GLuint id, Graph *graph) {
	glNewList(id, GL_COMPILE);
	
	glBegin(GL_POINTS);
	for (int i=0; i<(int)graph->nodes.size(); i++) {
		glVertex3dv(graph->nodes[i]->position);
	}
	glEnd();

	glEndList();
}











//////////////////////////////////////////////////////////////////////////////////////////////


//				BLUE NOISE, DO NOT TOUCH


//////////////////////////////////////////////////////////////////////////////////////////////









VolumeGrid::VolumeGrid(Volume *v, double ld, initial3d::random *r) {
	linkDistance = ld;
	double edgeSize = linkDistance / sqrt(3.0);
	baseSize = vec3d::one() * edgeSize;
	m_rand = r;
	volume = v;

	vec3d min = volume->getMin();
	vec3d max = volume->getMax();

	xSize = (max.x() - min.x()) / edgeSize;
	ySize = (max.y() - min.y()) / edgeSize;
	zSize = (max.z() - min.z()) / edgeSize;

	pointGrid = new Node *[xSize * ySize * zSize];
	cellGrid = new Cell *[xSize * ySize * zSize];
	indexOfVector = new int[xSize * ySize * zSize];

	// cout << "Creating from size : " << xSize << " : " << ySize << " : " << zSize << endl;

	for (int x = 0; x < xSize; x++) {
		for (int y = 0; y < ySize; y++) {
			for (int z = 0; z < zSize; z++) {
				int index = (z * ySize * xSize +  y * xSize + x);
				vec3d position(
					min.x()+(x * edgeSize), 
					min.y()+(y * edgeSize), 
					min.z()+(z * edgeSize));

				if (volume->contains(position)) {
					Cell *cell = new Cell(x, y, z, position, 0);
					active.insert(cell);
					cellGrid[index] = cell;
				} else {
					cellGrid[index] = NULL;
				}
				pointGrid[index] = NULL;
				indexOfVector[index] = -1;
			}
		}
	}
}

VolumeGrid::~VolumeGrid() {
	delete [] pointGrid;
	delete [] cellGrid;
	delete [] indexOfVector;

	for (int i=0; i<(int)children.size(); i++)
		delete children.at(i);
	for (set<Cell *>::iterator it = active.begin(); it != active.end(); ++it)
		delete *it;
}

void VolumeGrid::populate(Graph *graph) {

	double A = 0.7;
	int iteration = 0;

	while (!active.empty() > 0 && iteration < 1) {
		iteration++;
		// cout << "Iteration : " << iteration++ << endl;

		int estDarts = (int)(active.size() * A);
		set<Cell *> intersected;

		// cout << "estDarts" << estDarts << endl;
		// int i = 0;
		//throw darts
		for (int i = 0; i < estDarts && !active.empty(); i++) {
		// while (!active.empty()) {
			int index = m_rand->nextInt(active.size());
			// i++;

			Cell *c = NULL;
			set<Cell *>::iterator it = active.begin();
			for ( int l = 0; l < index && it != active.end(); ++it ) ++l;

			c = *it;

			double depthScale = 1.0/ math::pow(2.0, c->depth);
			vec3d minimum = c->position;
			vec3d scale = baseSize * depthScale;

			double x = m_rand->nextDouble() * scale.x() + minimum.x();
			double y = m_rand->nextDouble() * scale.y() + minimum.y();
			double z = m_rand->nextDouble() * scale.z() + minimum.z();
			vec3d point(x,y,z);

			// cout << "trying point " << i << " : " << point << endl;// << endl;

			if (!volume->contains(point)) {
				continue;
			}

			//if valid, add the nodes
			if (graph->validPoint(point, linkDistance)) { //TODO check graph 
				// cout << "valid point " << endl;
				// cout << "trying valid point " << i << " : " << point << endl;// << endl;

				//todo place node
				graph->addAndLinkNode(new Node(point), linkDistance * 3); //TODO check graph 
				//remove this node as it cannot be chosen again
				active.erase(c);

				updateActive(point, c, intersected);

				//remove all of the active nodes the are encompassed
				//and add the nodes that are intersected to
				//THIS REALLY NEEDS TO BE OPTIMISED TODO
				// for (vector<Cell *>::iterator it = active.begin(); it != active.end(); ) {
				// 	//TODO check intersect and flag
				// 	if (cellInside(*it, point))
				// 		it = active.erase(it);
				// 	else if (cellInside(*it)) {
				// 		intersected.insert(*it);
				// 		it = active.erase(it);
				// 	} else
				// 		++it;
				// }
			}
			// cout << "end " << endl;
		}
		// set<Cell *> tempActive;

		// cout << "Iterate " << endl;

		// //iterate
		// for (set<Cell *>::iterator it = active.begin(); it != active.end(); ++it) {
		// 	// cout << "Spliting cell " << endl;
		// 	splitCell(*it, tempActive);
		// }
		// for (set<Cell *>::iterator it = intersected.begin(); it != intersected.end(); ++it) {
		// 	// cout << "Spliting and check cell " << endl;
		// 	splitAndCheckCell(*it, tempActive);
		// }
		// active = tempActive;

		// cout << "asdasdasdasdasd" << endl;
	}
}




void VolumeGrid::updateActive(vec3d point, Cell *cell, set<Cell *> &intersected) {
	for (int x = std::max(cell->x-1, 0); x<std::min(xSize, cell->x+1); x++) {
		for (int y = std::max(cell->y-1, 0); y<std::min(ySize, cell->y+1); y++) {
			for (int z = std::max(cell->z-1, 0); z<std::min(zSize, cell->z+1); z++) {

				int index = (z * ySize * xSize +  y * xSize + x);
				int vIndex = indexOfVector[index];

				if (children.size() == 0) {
					Cell * c = cellGrid[index];
					if (c != NULL) {
						// if (cellIntersects(c, point)) {
							// active.erase(c);
							// intersected.insert(c);
						// } else 
						if (cellInside(c, point)) {
							active.erase(c);
						}
					}
				} else if (vIndex != -1) {
					set<Cell *>::iterator it = children.at(vIndex)->begin();
					for (; it != children.at(vIndex)->end(); ++it) {
						// if (cellIntersects(*it, point)) {
							// active.erase(*it);
							// intersected.insert(*it);
						// } else 
						if (cellInside(*it, point)) {
							active.erase(*it);
						}
					}
				}
			}
		}
	}
}


//split and add children to the active list
void VolumeGrid::splitCell(Cell *c, set<Cell *> &active) {
	double depthScale = 1.0/ math::pow(2.0, c->depth);
	vec3d minimum = c->position;
	vec3d scale = (baseSize * depthScale) * 0.5;
	
	double x = minimum.x();
	double y = minimum.y();
	double z = minimum.z();

	double sx = scale.x();
	double sy = scale.y();
	double sz = scale.z();

	Cell * cell;

	int index = (z * ySize * xSize +  y * xSize + x);
	int vIndex = indexOfVector[index];
	set<Cell *> *child;

	if (vIndex < 0) {
		indexOfVector[index] = children.size();
		child = new set<Cell *>;
		children.push_back(child);
	} else {
		child = children.at(vIndex);
		child->erase(c);
	}

	cell = new Cell(c->x, c->y, c->z, vec3d(   x,    y,    z), c->depth+1);
	active.insert(cell); child->insert(cell);
	cell = new Cell(c->x, c->y, c->z, vec3d(x+sx,    y,    z), c->depth+1);
	active.insert(cell); child->insert(cell);
	cell = new Cell(c->x, c->y, c->z, vec3d(   x, y+sy,    z), c->depth+1);
	active.insert(cell); child->insert(cell);
	cell = new Cell(c->x, c->y, c->z, vec3d(   x,    y, z+sz), c->depth+1);
	active.insert(cell); child->insert(cell);
	cell = new Cell(c->x, c->y, c->z, vec3d(x+sx, y+sy,    z), c->depth+1);
	active.insert(cell); child->insert(cell);
	cell = new Cell(c->x, c->y, c->z, vec3d(x+sx,    y, z+sz), c->depth+1);
	active.insert(cell); child->insert(cell);
	cell = new Cell(c->x, c->y, c->z, vec3d(   x, y+sy, z+sz), c->depth+1);
	active.insert(cell); child->insert(cell);
	cell = new Cell(c->x, c->y, c->z, vec3d(x+sx, y+sy, z+sz), c->depth+1);
	active.insert(cell); child->insert(cell);

	delete c;
}

//split and add children to the active list
//delete pointer when finished
void VolumeGrid::splitAndCheckCell(Cell *c, set<Cell *> &active) {

	double depthScale = 1.0/ math::pow(2.0, c->depth);
	vec3d minimum = c->position;
	vec3d scale = (baseSize * depthScale) * 0.5;
	
	double x = minimum.x();
	double y = minimum.y();
	double z = minimum.z();

	double sx = scale.x();
	double sy = scale.y();
	double sz = scale.z();

	Cell * cell;

	int index = (z * ySize * xSize +  y * xSize + x);
	int vIndex = indexOfVector[index];
	set<Cell *> *child;

	if (vIndex >= children.size())
		vIndex = -1;

	if (vIndex < 0) {
		// cout << " No Exists a stuff " << children.size() << " index  " << vIndex << endl;
		indexOfVector[index] = children.size();
		child = new set<Cell *>;
		children.push_back(child);

	} else {
		// cout << "Exists a stuff " << children.size() << " index  " << vIndex << endl;
		child = children.at(vIndex);
		child->erase(c);
	}

	if (!cellInside( cell = new Cell(c->x, c->y, c->z, vec3d(   x,    y,    z), c->depth+1) )) {
		active.insert(cell); child->insert(cell); }
	if (!cellInside( cell = new Cell(c->x, c->y, c->z, vec3d(x+sx,    y,    z), c->depth+1) )) {
		active.insert(cell); child->insert(cell); }
	if (!cellInside( cell = new Cell(c->x, c->y, c->z, vec3d(   x, y+sy,    z), c->depth+1) )) {
		active.insert(cell); child->insert(cell); }
	if (!cellInside( cell = new Cell(c->x, c->y, c->z, vec3d(   x,    y, z+sz), c->depth+1) )) {
		active.insert(cell); child->insert(cell); }
	if (!cellInside( cell = new Cell(c->x, c->y, c->z, vec3d(x+sx, y+sy,    z), c->depth+1) )) {
		active.insert(cell); child->insert(cell); }
	if (!cellInside( cell = new Cell(c->x, c->y, c->z, vec3d(x+sx,    y, z+sz), c->depth+1) )) {
		active.insert(cell); child->insert(cell); }
	if (!cellInside( cell = new Cell(c->x, c->y, c->z, vec3d(   x, y+sy, z+sz), c->depth+1) )) {
		active.insert(cell); child->insert(cell); }
	if (!cellInside( cell = new Cell(c->x, c->y, c->z, vec3d(x+sx, y+sy, z+sz), c->depth+1) )) {
		active.insert(cell); child->insert(cell); }

	delete c;
}







bool VolumeGrid::validPoint(vec3d point, Cell *cell) {
	for (int x = std::max(cell->x-1, 0); x<std::min(xSize, cell->x+1); x++) {
		for (int y = std::max(cell->y-1, 0); y<std::min(ySize, cell->y+1); y++) {
			for (int z = std::max(cell->z-1, 0); z<std::min(zSize, cell->z+1); z++) {

				int index = (z * ySize * xSize +  y * xSize + x);
				Node *otherPoint = pointGrid[index];

				if (otherPoint != NULL) {
					if (+(otherPoint->position - point) < linkDistance) {
						return false;
					}
				}
			}
		}
	}
	return true;
}



bool VolumeGrid::cellInside(Cell *cell, vec3d point) {
	double resultX = 0, resultY = 0, resultZ = 0;

	double scale = 1.0/ math::pow(2.0, cell->depth);
	vec3d minimum = cell->position;
	vec3d maximum = cell->position + (baseSize * scale);

	vec3d center = point;

	resultX = std::max(maximum.x() - center.x(), center.x() - minimum.x());
	resultY = std::max(maximum.y() - center.y(), center.y() - minimum.y());
	resultZ = std::max(maximum.z() - center.z(), center.z() - minimum.z());

	vec3d farestPoint(resultX, resultY, resultZ);
	return linkDistance >= +farestPoint;
}


bool VolumeGrid::cellInside(Cell *cell) {
	double resultX = 0, resultY = 0, resultZ = 0;

	double scale = 1.0/ math::pow(2.0, cell->depth);
	vec3d minimum = cell->position;
	vec3d maximum = cell->position + (baseSize * scale);

	for (int x = std::max(cell->x-1, 0); x<std::min(xSize, cell->x+1); x++){
		for (int y = std::max(cell->y-1, 0); y<std::min(ySize, cell->y+1); y++){
			for (int z = std::max(cell->z-1, 0); z<std::min(zSize, cell->z+1); z++){

				int index = (z * ySize * xSize +  y * xSize + x);
				Node *otherPoint = pointGrid[index];

				if (otherPoint != NULL) {

					double resultX = 0, resultY = 0, resultZ = 0;
					vec3d center = otherPoint->position;

					resultX = std::max(maximum.x() - center.x(), center.x() - minimum.x());
					resultY = std::max(maximum.y() - center.y(), center.y() - minimum.y());
					resultZ = std::max(maximum.z() - center.z(), center.z() - minimum.z());

					vec3d farestPoint(resultX, resultY, resultZ);
					if (linkDistance >= +farestPoint)
					return true;
				}
			}
		}
	}
	return false;
}



bool VolumeGrid::cellIntersects(Cell *cell, vec3d point) {
	using namespace initial3d;

	double scale = 1.0/ math::pow(2.0, cell->depth);
	vec3d minimum = cell->position;
	vec3d maximum = cell->position + (baseSize * scale);

	double resultX = 0, resultY = 0, resultZ = 0;
	vec3d center = point;

	// for each plane record the 1-d vector from the
	// edges to the sphere on that plane.
	if (center.x() < minimum.x())
		resultX = center.x() - minimum.x();
	else if (center.x() > maximum.x())
		resultX = center.x() - maximum.x();

	if (center.y() < minimum.y())
		resultY = center.y() - minimum.y();
	else if (center.y() > maximum.y())
		resultY = center.y() - maximum.y();

	if (center.z() < minimum.z())
		resultY = center.z() - minimum.z();
	else if (center.z() > maximum.z())
		resultY = center.z() - maximum.z();

	// compile the vector together and check the distance
	vec3d collisionNorm(resultX, resultY, resultZ);
	return linkDistance >= +collisionNorm;
}


// bool VolumeGrid::cellIntersects(Cell *cell) {
// 	using namespace initial3d;

// 	double scale = 1.0/ math::pow(2, cell->depth);
// 	vec3d minimum = cell->position;
// 	vec3d maximum = cell->position + (baseSize * scale);


// 	for (int x = std::max(cell->x-2, 0); x<std::min(xSize, cell->x+2); x++){
// 		for (int y = std::max(cell->y-2, 0); y<std::min(ySize, cell->y+2); y++){
// 			for (int z = std::max(cell->z-2, 0); z<std::min(zSize, cell->z+2); z++){

// 				int index = (z * ySize * xSize +  y * xSize + x);
// 				Cell *nCell = cellGrid[index];

// 				if (nCell != NULL && nCell.node != NULL) {

// 					double resultX = 0, resultY = 0, resultZ = 0;
// 					Node *node = nCell->node;
// 					vec3d center = node->position;

// 					// for each plane record the 1-d vector from the
// 					// edges to the sphere on that plane.
// 					if (center.x() < minimum.x())
// 						resultX = center.x() - minimum.x();
// 					else if (center.x() > maximum.x())
// 						resultX = center.x() - maximum.x();

// 					if (center.y() < minimum.y())
// 						resultY = center.y() - minimum.y();
// 					else if (center.y() > maximum.y())
// 						resultY = center.y() - maximum.y();

// 					if (center.z() < minimum.z())
// 						resultY = center.z() - minimum.z();
// 					else if (center.z() > maximum.z())
// 						resultY = center.z() - maximum.z();

// 					// compile the vector together and check the distance
// 					vec3d collisionNorm(resultX, resultY, resultZ);
// 					if (linkDistance >= +collisionNorm)
// 						return true;
// 				}
// 			}
// 		}
// 	}
// 	return false;
// }
