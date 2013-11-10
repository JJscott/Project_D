
#pragma once

#include <map>
#include <iostream>
#include <vector>
#include <set>
#include "GLee.h"
#include <GL/glu.h>
#include "volume.h"
#include "initial3d.h"
#include "octree.h"


struct Edge;
struct Node {
	initial3d::vec3d position;
	std::vector<Edge *> edges;
	typedef double coord_t;

	initial3d::aabbd getAABB () const {
		return initial3d::aabbd(position, position);
	}

	Node(initial3d::vec3d p) : position(p) {}
};
struct NodeP{
	Node *node;
	double cost;
	bool operator<(const NodeP &rhs) const {
		return cost > rhs.cost;
	}
	NodeP(Node *n, double cost_) : node(n), cost(cost_) {}
};

struct Edge {
	Node *n1, *n2;
	double weight;
	Node * getOther(Node *n) { return n==n2 ? n1 : n2; }
	Edge(Node *nn1, Node *nn2) : n1(nn1), n2(nn2), weight(0) {
		n1->edges.push_back(this);
		n2->edges.push_back(this);
	}
	Edge(Node *nn1, Node *nn2, double w) : n1(nn1), n2(nn2), weight(w) {
		n1->edges.push_back(this);
		n2->edges.push_back(this);
	}
};






struct Graph {
	std::vector<Node *> nodes;
	std::vector<Edge *> edges;
	initial3d::octree<Node> octopuss;

	bool validPoint(initial3d::vec3d point, double minDistance){
		std::vector<Node *> found;
		octopuss.find(found, initial3d::aabbd::fromOuterSphere(point, minDistance));
		return found.empty();
	}

	void addAndLinkNode(Node *node, double linkDistance) {

		std::vector<Node *> found;
		octopuss.find(found, initial3d::aabbd::fromInnerSphere(node->position, linkDistance));

		// Hi josh, how are you?                                                                                                                                                                todo check distance proparly

		for (int i=0; i<(int)found.size(); i++) {
			Edge * edge = new Edge(node, found[i], +(node->position - found[i]->position));
			edges.push_back(edge);
			node->edges.push_back(edge);
			found[i]->edges.push_back(edge);
		}



		octopuss.add(node);
		nodes.push_back(node);
	}


	~Graph() {
		for (int i=0; i<(int)nodes.size(); i++) {
			delete nodes[i];
		}
		for (int i=0; i<(int)edges.size(); i++) {
			delete edges[i];
		}
	}
};




struct TreeNode {
	int branchCount;
	double radius;
	initial3d::vec3d position;
	std::vector<TreeNode *> children;

	TreeNode(initial3d::vec3d p) : branchCount(0), radius(0), position(p) {}

	int calculateBranchCount() {
		int count = 1;
		for (int i = 0; i < (int)children.size(); i++) {
			count += children[i]->calculateBranchCount();
		}
		branchCount = count;
		return count;
	}

	void calculateRadius(double newRadius, double exponent) {
		using namespace initial3d;
		radius = newRadius;

		double area = math::pi() * math::pow(radius, exponent);

		for (int i = 0; i < (int)children.size(); i++) {
			double ratio = children.at(i)->branchCount / (double)branchCount;
			double childradius = math::sqrt((ratio * area) / math::pi()); //todo change sqrt to exponent root
			children.at(i)->calculateRadius(childradius, exponent);
		}
	}
};

struct TempTreeNode{
	Node *node;
	std::vector<TempTreeNode*> children;

	TreeNode * toTreeNode() {
		TreeNode *tree = new TreeNode(node->position);
		for (int i=0; i<(int)children.size(); i++) {
			tree->children.push_back(children.at(i)->toTreeNode());
		}
		return tree;
	}

	TempTreeNode(Node *n) : node(n) {}
};







class TreeGenerator{
private:

	double param_nodeLinkDistance;
	double param_subgraphSize;
	double param_subgraphAngle;

	int m_maxDepth;
	int m_branchFactor; //number of endpoints per volume
	int m_graphPopulation; //number of nodes in each graph


	initial3d::random *m_rand;

public:
	TreeGenerator(); //temp
	// TreeGenerator(graphResolution, linkingDistance, shrinkingParameter, sizingAngle);
	void buildTree(GLuint);
	void buildSubTree(TempTreeNode *, std::vector<TempTreeNode *> &, int);

	void createGraph(vol::Volume *, Graph *);
	void createGraph2(vol::Volume *, Graph *);
	void mapParents(Graph *, Node *, std::map<Node *, Node *> &);

	// TreeNode * buildTree(Graph *, Node *, std::map<Node *, Node *> &, std::vector<Node*> &);
	void buildTree(Graph *, TempTreeNode *, std::map<Node *, Node *> &, std::vector<Node*> &, std::vector<TempTreeNode *> &);
	TreeNode * buildNaiveTree(Graph *, Node *, std::map<Node *, Node *> &);
	void createMesh(GLuint, TreeNode *);
	void drawBranch(GLUquadric *, TreeNode *);
	void drawCylinder(GLUquadric *, const initial3d::vec3d &, double, const initial3d::vec3d &, double);


	void drawGraphPoints(GLuint, Graph *);
};


struct Cell {
	int x, y, z;
	initial3d::vec3d position;
	int depth;

	Cell(int xs, int ys, int zs, initial3d::vec3d pos, int d) : x(xs), y(ys), z(zs), position(pos), depth(d) {}
};


class VolumeGrid {
private:
	int xSize;
	int ySize;
	int zSize;

	Node **pointGrid;
	Cell **cellGrid;
	int *indexOfVector;
	std::vector<std::set<Cell *> *> children;

	vol::Volume *volume;

	initial3d::vec3d baseSize;
	std::set<Cell *> active;

	double linkDistance;
	initial3d::random* m_rand;


	void updateActive(initial3d::vec3d point, Cell *cell, std::set<Cell *> &intersected);
	void splitCell(Cell *c, std::set<Cell *> &active);
	void splitAndCheckCell(Cell *c, std::set<Cell *> &active);
	bool validPoint(initial3d::vec3d point,  Cell *cell);
	bool cellInside(Cell *cell, initial3d::vec3d point);
	bool cellInside(Cell *cell);
	bool cellIntersects(Cell *cell, initial3d::vec3d point);

public:

	VolumeGrid(vol::Volume *volume, double ld, initial3d::random *r);
	~VolumeGrid();
	void populate(Graph *graph);
};