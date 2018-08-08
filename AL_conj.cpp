#include<iostream>
#include<vector>
#include<algorithm>
#include<fstream>
#include<tuple>

extern "C" {
	#include "glpk-4.60/src/glpk.h"
}

std::ofstream outfile("graphs/list.txt");

constexpr int DEGREE = 3;
constexpr int MAX_RADIUS = 2;
constexpr int COLORS = 2;
constexpr double PROB = 1.0/DEGREE;
constexpr double EPSILON = 0.000001;

constexpr int MAX_CONSTRAINTS = 1000002;

// coefficients for linear programming matrix
int ia[MAX_CONSTRAINTS], ja[MAX_CONSTRAINTS];
double ar[MAX_CONSTRAINTS];

// colored trees will be represented as a tuple of indices of trees of smaller radius
struct Tree {
	int radius;
	int root_color;
	std::vector<int> children_types;

	// I also want to remember what tree is this if I cut off the leaves ...
	// this value is index in (ir)regular_trees[radius - 1]
	int cut_off_type;

	const Tree& get_subtree(int k);

	Tree(int r, int r_c, int c_o_t, std::vector<int> ch_t) :
	  children_types(std::move(ch_t)),
		radius(r),
		root_color(r_c),
		cut_off_type(c_o_t) {}

	bool operator==(const Tree& lhs) const {
		return lhs.radius == radius &&
		       lhs.root_color == root_color &&
					 lhs.children_types == children_types;
	}
};

using EdgeTree = std::pair<Tree, Tree>;
using PathTree = std::tuple<Tree, Tree, Tree>;

// trees[k] are trees with root of degree k, the rest of the vertices have degree DEGREE
std::vector<Tree> trees[DEGREE + 1][MAX_RADIUS + 1];
std::vector<EdgeTree> e_trees[MAX_RADIUS + 1];
std::vector<PathTree> p_trees[MAX_RADIUS + 1];

const Tree& Tree::get_subtree(int k) {
	// This will only be used with regular trees
	if((children_types.size() != DEGREE) || k < 0 || k > DEGREE)
		std::cerr << "Tree get_subtree error.\n";
	return trees[DEGREE - 1][radius - 1][children_types[k]];
}

void write_tree(const Tree& t) {
	std::vector<Tree> v1, v2;
	int r = t.radius;
	v2.push_back(std::move(t));
	do {
		v1 = std::move(v2);
		v2.clear();
		outfile << "  ";

		for(size_t j = 0; j < v1.size(); j++) {
			outfile << v1[j].root_color;
			for(size_t k = 0; k < v1[j].children_types.size(); k++) {
				v2.push_back(trees[DEGREE - 1][r - 1][v1[j].children_types[k]]);
			}
		}
		r--;
	} while(!v2.empty());

	outfile << ",  cut off type: " << t.cut_off_type << '\n';
}

void write_to_file() {
	for(int r = 0; r <= MAX_RADIUS; r++) {
		for(int deg = 1; deg <= DEGREE; deg++) {
			outfile << "\nVERTEX-ROOTED TREES, RADIUS " << r <<
			           ", ROOT DEGREE " << deg << ":\n";
			for(size_t i = 0; i < trees[deg][r].size(); i++) {
				outfile << i + 1 << ":";
				write_tree(trees[deg][r][i]);
			}
		}
	}

	outfile << "\nREGULAR EDGE-ROOTED TREES, RADIUS " << MAX_RADIUS - 1 << ":\n";
	for(size_t i = 0; i < e_trees[MAX_RADIUS - 1].size(); i++) {
		outfile << i + 1 << ":\t";
		write_tree(e_trees[MAX_RADIUS - 1][i].first);
		outfile << "\t";
		write_tree(e_trees[MAX_RADIUS - 1][i].second);
	}

	outfile << "\nREGULAR PATH-ROOTED TREES, RADIUS " << MAX_RADIUS - 1 << ":\n";
	for(size_t i = 0; i < p_trees[MAX_RADIUS - 1].size(); i++) {
		outfile << i + 1 << ":\t";
		write_tree(std::get<0>(p_trees[MAX_RADIUS - 1][i]));
		outfile << "\t";
		write_tree(std::get<1>(p_trees[MAX_RADIUS - 1][i]));
		outfile << "\t";
		write_tree(std::get<2>(p_trees[MAX_RADIUS - 1][i]));
	}
}

// find the tree in the vector trees with desired children types
// return -1 if not found
int find_tree(std::vector<int>&& children_types,
              int root_color, int degree, int radius) {
	std::sort(children_types.begin(), children_types.end());
	for(size_t i = 0; i < trees[degree][radius].size(); i++) {
		if(trees[degree][radius][i].children_types == children_types &&
		    trees[degree][radius][i].root_color == root_color)
			return i;
	}
	std::cerr << "Find vertex-rooted tree error. Degree: " << degree <<
	             ", Radius: " << radius << "\n";
	exit(1);
	return -1;
}

// the rooting directed edge is from t1 to t2
int find_edge_tree(int a, int b) {
	for(size_t i = 0; i < e_trees[MAX_RADIUS - 1].size(); i++) {
		if(e_trees[MAX_RADIUS - 1][i].first == trees[DEGREE - 1][MAX_RADIUS - 1][a] &&
				e_trees[MAX_RADIUS - 1][i].second == trees[DEGREE - 1][MAX_RADIUS - 1][b])
			return i;
	}
	std::cerr << "Find edge-rooted tree error.\n";
	exit(1);
	return -1;
}

int find_path_tree(int a, int b, int c) {
	if(a >= trees[DEGREE - 1][MAX_RADIUS - 2].size() ||
	    c >= trees[DEGREE - 1][MAX_RADIUS - 2].size() ||
			b >= trees[DEGREE - 2][MAX_RADIUS - 1].size()) {
		std::cerr << "Find path-rooted tree parameter error.\n";
		exit(1);
	}
	for(size_t i = 0; i < p_trees[MAX_RADIUS - 1].size(); i++) {
		if(std::get<0>(p_trees[MAX_RADIUS - 1][i]) == trees[DEGREE - 1][MAX_RADIUS - 2][a] &&
				std::get<1>(p_trees[MAX_RADIUS - 1][i]) == trees[DEGREE - 2][MAX_RADIUS - 1][b] &&
				std::get<2>(p_trees[MAX_RADIUS - 1][i]) == trees[DEGREE - 1][MAX_RADIUS - 2][c])
			return i;
	}
	std::cerr << "Find path-rooted tree error.\n";
	exit(1);
	return -1;
}

// parameter is a regular tree and one of the neighbours of its root, return index in irregular trees after cutting away the whole subtree, final parameter is where to look for it
int kill_branch(const Tree& t, int nb_type, int degree, int radius) {
	std::vector<int> vec = t.children_types;
	auto p = std::find(vec.begin(), vec.end(), nb_type);
	if(p == vec.end()) { std::cerr << "Kill branch error.\n"; exit(1); }
	vec.erase(p);
	return find_tree(std::move(vec), t.root_color, degree, radius);
}

// helper container for function recurse
std::vector<int> children;
void recurse(size_t p, std::vector<Tree>& vec, int radius, int root_degree) {
	// we can assume radius >= 1
	if(p == root_degree) {
		for(int i = 0; i < COLORS; i++) {
			// i is the color of root
			// I need to calculate the cut_off_type
			std::vector<int> temp;
			for(size_t j = 0; j < children.size(); j++) {
				// find the right cut off type of this tree of radius = radius - 1
				temp.push_back(trees[DEGREE - 1][radius - 1][children[j]].cut_off_type);
			}
			int cut_off_type;
			if(radius == 1) cut_off_type = i;
			else {
				cut_off_type = find_tree(std::move(temp), i, root_degree, radius - 1);
			}
			Tree t(radius, i, cut_off_type, children);
			vec.push_back(t);
		}
	} else {
		int start = p == 0 ? 0 : children[p-1];
		for(size_t i = start; i < trees[DEGREE - 1][radius - 1].size(); i++) {
			children[p] = i;
			recurse(p + 1, vec, radius, root_degree);
		}
	}
}

void generate_graphs(int r) {
	if(r > MAX_RADIUS) return;
	std::vector<Tree> vec;

	if(r == 0) {
		for(int i = 0; i < COLORS; i++) {
			// just root with color c should have index c
			vec.emplace_back(0, i, -1, std::vector<int>());
		}
		for(int deg = 1; deg <= DEGREE; deg++) trees[deg][r] = vec;
	} else {
		// I dont need to generate irregular graphs  graphs of maximal radius:
		for(int deg = 1; deg <= DEGREE; deg++) {
			children.clear();
			children.resize(deg);
			vec.clear();
			recurse(0, vec, r, deg);
			trees[deg][r] = std::move(vec);
			std::cout << "Generated " << trees[deg][r].size() <<
				" regular vertex-rooted trees with root-degree " << deg <<
				" and radius " << r << ".\n";
		}
	}

	if(r == MAX_RADIUS) {
		// edge-rooted trees:
		for(size_t i = 0; i < trees[DEGREE - 1][MAX_RADIUS - 1].size(); i++)
		for(size_t j = 0; j < trees[DEGREE - 1][MAX_RADIUS - 1].size(); j++) {
			e_trees[MAX_RADIUS - 1].emplace_back(
				trees[DEGREE - 1][MAX_RADIUS - 1][i], trees[DEGREE - 1][MAX_RADIUS - 1][j]
			);
		}
		std::cout << "Generated " << e_trees[MAX_RADIUS - 1].size() <<
		             " regular edge-rooted trees of radius " << MAX_RADIUS - 1 << ".\n";

		// path-rooted trees:
		for(size_t i = 0; i < trees[DEGREE - 1][MAX_RADIUS - 2].size(); i++)
		for(size_t j = 0; j < trees[DEGREE - 2][MAX_RADIUS - 1].size(); j++)
		for(size_t k = 0; k < trees[DEGREE - 1][MAX_RADIUS - 2].size(); k++) {
			p_trees[MAX_RADIUS - 1].emplace_back(
				trees[DEGREE - 1][MAX_RADIUS - 2][i],
				trees[DEGREE - 2][MAX_RADIUS - 1][j],
				trees[DEGREE - 1][MAX_RADIUS - 2][k]);
	 	}
		std::cout << "Generated " << p_trees[MAX_RADIUS - 1].size() <<
		             " regular path-rooted trees of radius " << MAX_RADIUS - 1 << ".\n";
	}

	generate_graphs(r + 1);
}

// returns true if d is very very small
bool check_epsilon(double d) {
	return d <= EPSILON && d >= (-EPSILON);
}

void check_epsilon_ar(int length) {
	int counter = 0;
	for(int i = 1; i <= length; i++) {
		if(check_epsilon(ar[i])) {
			ar[i] = 0;
			counter++;
		}
	}
	std::cout << "Removed " << counter << " epsilons.\n";
}

void write_to_matrix(int row_index, int column_index, double val, int& counter) {
	ia[counter] = row_index;
	ja[counter] = column_index;
	ar[counter] = val;
	counter++;
}

void sum_to_one(int& counter, int row_index, int column_count) {
	for(int i = 0; i < column_count; i++)
		write_to_matrix(row_index, i + 1, 1.0, counter);
}

int generate_constraints_rerooting() {
	int counter = 1;
	int start_counter;
	int * p;
	int row_index;
	int temp;
	int new_root_color;
	for(size_t i = 0; i < trees[DEGREE][MAX_RADIUS].size(); i++) {
		// first, consider just cutting off the leaves
		// INDEXING IS SHIFTED!!!
		start_counter = counter;
		//std::cout << "start_counter: " << start_counter << "  value of i: " << i <<'\n';
		ia[counter] = trees[DEGREE][MAX_RADIUS][i].cut_off_type + 1;
		ja[counter] = i + 1;
		ar[counter] = (-1);
		counter++;

		// second, reroot to random neighbour (the graphs is regular, so we reroot uniformly)
		for(int nb = 0; nb < DEGREE; nb++) {
			// find color of the new root
			new_root_color = trees[DEGREE][MAX_RADIUS][i].get_subtree(nb).root_color;
			// get children of the subtree:
			std::vector<int> vec(trees[DEGREE][MAX_RADIUS][i].get_subtree(nb).children_types);
			// now cut off leaves
			const Tree& t_cut =
				trees[DEGREE][MAX_RADIUS - 1][trees[DEGREE][MAX_RADIUS][i].cut_off_type];
			// now cut off the branch
			temp = trees[DEGREE][MAX_RADIUS][i].children_types[nb];
			temp = kill_branch(t_cut,
				trees[DEGREE - 1][MAX_RADIUS - 1][temp].cut_off_type,
				DEGREE - 1, MAX_RADIUS - 1);
			// finally cut again
			temp = trees[DEGREE - 1][MAX_RADIUS - 1][temp].cut_off_type;
			vec.push_back(temp);

			row_index = find_tree(std::move(vec), new_root_color, DEGREE, MAX_RADIUS - 1) + 1;
			// check if this has value has been set before
			// some pointer arithmetic magic
			p = std::find(ia + start_counter, ia + counter, row_index);
			if(p == ia + counter) {
				write_to_matrix(row_index, i + 1, PROB, counter);
			} else {
				// more pointer arithmetic magic
				if(ja[p - ia] != i + 1) {
					std::cout << p - ia << '\n';
					std::cout << "Error arithmetic magic! Rerooting. ja value: " <<
						ja[p - ia - start_counter + 1] <<
						"   counter: " << counter << "   i value: " << i <<
						"\n"; exit(1);
				}
				ar[p - ia] += PROB;
			}
		}
	}
	check_epsilon_ar(counter - 1);
	sum_to_one(counter, trees[DEGREE][MAX_RADIUS - 1].size() + 1,
		trees[DEGREE][MAX_RADIUS].size());
	counter--;
	std::cout << "Rerooting polyhedron: Generated " << counter << " coefficients.\n";
	return counter;
}

int generate_constraints_AL() {
	int counter = 1;
	int start_counter;
	int temp;
	int subtree_index_1;
	int subtree_index_2;
	int row_index;
	int *p;
	for(size_t i = 0; i < trees[DEGREE][MAX_RADIUS].size(); i++) {
		start_counter = counter;
		for(int nb = 0; nb < DEGREE; nb++) {
			// find the subtree of chosen neighbour
			subtree_index_1 = trees[DEGREE][MAX_RADIUS][i].children_types[nb];

			subtree_index_2 = kill_branch(trees[DEGREE][MAX_RADIUS][i],
				trees[DEGREE][MAX_RADIUS][i].children_types[nb], DEGREE - 1, MAX_RADIUS);
			subtree_index_2 = trees[DEGREE - 1][MAX_RADIUS][subtree_index_2].cut_off_type;

			// one direction of edge:
			row_index = find_edge_tree(subtree_index_1, subtree_index_2) + 1;
			// check if this has value has been set before
			// some pointer arithmetic magic
			p = std::find(ia + start_counter, ia + counter, row_index);
			if(p == ia + counter) {
				write_to_matrix(row_index, i + 1, PROB, counter);
			} else {
				// more pointer arithmetic magic
				if(ja[p - ia] != i + 1) { std::cerr << "Error arithmetic magic!\n"; exit(1); }
				ar[p - ia] += PROB;
			}
			// the other direction of edge:
			row_index = find_edge_tree(subtree_index_2, subtree_index_1) + 1;
			p = std::find(ia + start_counter, ia + counter, row_index);
			if(p == ia + counter) {
				write_to_matrix(row_index, i + 1, (-1)*PROB, counter);
			} else {
				// more pointer arithmetic magic
				if(ja[p - ia] != i + 1) { std::cerr << "Error arithmetic magic!\n"; exit(1); }
				ar[p - ia] -= PROB;

			}

		}
	}
	check_epsilon_ar(counter - 1);
	sum_to_one(counter, e_trees[MAX_RADIUS - 1].size() + 1,
		trees[DEGREE][MAX_RADIUS].size());
	counter--;
	std::cout << "A-L polyhedron: Generated " << counter << " coefficients.\n";
	return counter;

}

int generate_constraints_path() {
	int counter = 1;
	int start_counter;
	int temp;
	int subtree_index_1, subtree_index_2, subtree_index_3;
	int row_index;
	int *p;
	for(size_t i = 0; i < trees[DEGREE][MAX_RADIUS].size(); i++) {
		start_counter = counter;
		for(int nb1 = 0; nb1 < DEGREE; nb1++) for(int nb2 = 0; nb2 < DEGREE - 1; nb2++) {
			//std::cout << counter << '\n';
			subtree_index_1 = trees[DEGREE][MAX_RADIUS][i].children_types[nb1];
			subtree_index_1 = trees[DEGREE - 1][MAX_RADIUS - 1][subtree_index_1].children_types[nb2];

			subtree_index_2 = trees[DEGREE][MAX_RADIUS][i].children_types[nb1];
			// THE FOLLOWING LINES SHOULD BE CORRECT
			subtree_index_2 =
				kill_branch(trees[DEGREE - 1][MAX_RADIUS - 1][subtree_index_2],
					subtree_index_1, DEGREE - 2, MAX_RADIUS - 1);

			subtree_index_3 =
				kill_branch(trees[DEGREE][MAX_RADIUS][i],
					trees[DEGREE][MAX_RADIUS][i].children_types[nb1], DEGREE - 1, MAX_RADIUS);
			subtree_index_3 = trees[DEGREE - 1][MAX_RADIUS][subtree_index_3].cut_off_type;
			subtree_index_3 = trees[DEGREE - 1][MAX_RADIUS - 1][subtree_index_3].cut_off_type;

			// DIRECTION OF EDGES OUT OF ROOT:
			row_index = find_path_tree(subtree_index_1, subtree_index_2, subtree_index_3);
			// check if this has value has been set before
			// some pointer arithmetic magic
			p = std::find(ia + start_counter, ia + counter, 3 * row_index + 1);
			if(p == ia + counter) {
				write_to_matrix(3 * row_index + 1, i + 1, PROB, counter);
			} else {
				// more pointer arithmetic magic
				if(ja[p - ia] != i + 1) { std::cerr << "Error arithmetic magic! Path11\n"; exit(1); }
				ar[p - ia] += PROB;
			}
			p = std::find(ia + start_counter, ia + counter, 3 * row_index + 2);
			if(p == ia + counter) {
				write_to_matrix(3 * row_index + 2, i + 1, PROB, counter);
		 	} else {
				if(ja[p - ia] != i + 1) { std::cerr << "Error arithmetic magic! Path12\n"; exit(1); }
				ar[p - ia] += PROB;
			}

			// DIRECTION OF EDGE TOWARDS ROOT:
			row_index = find_path_tree(subtree_index_3, subtree_index_2, subtree_index_1);

			p = std::find(ia + start_counter, ia + counter, 3 * row_index + 1);
			if(p == ia + counter) {
				write_to_matrix(3 * row_index + 1, i + 1, (-1)*PROB, counter);
			} else {
				// more pointer arithmetic magic
				if(ja[p - ia] != i + 1) { std::cerr << "Error arithmetic magic! Path21\n"; exit(1); }
				ar[p - ia] -= PROB;
			}
			p = std::find(ia + start_counter, ia + counter, 3 * row_index + 3);
			if(p == ia + counter) {
				write_to_matrix(3 * row_index + 3, i + 1, (-1)*PROB, counter);
		 	} else {
				if(ja[p - ia] != i + 1) { std::cerr << "Error arithmetic magic! Path22\n"; exit(1); }
				ar[p - ia] -= PROB;
			}

			// THIRD POSSIBILITY:
			subtree_index_1 = trees[DEGREE][MAX_RADIUS][i].children_types[nb1];
			subtree_index_1 = trees[DEGREE - 1][MAX_RADIUS - 1][subtree_index_1].cut_off_type;

			subtree_index_2 = trees[DEGREE][MAX_RADIUS][i].cut_off_type;
			subtree_index_2 =
				kill_branch(trees[DEGREE][MAX_RADIUS - 1][subtree_index_2],
					trees[DEGREE - 1][MAX_RADIUS - 1][trees[DEGREE][MAX_RADIUS][i].children_types[nb1]].cut_off_type,
					DEGREE - 1, MAX_RADIUS - 1);
			subtree_index_2 =
				kill_branch(trees[DEGREE - 1][MAX_RADIUS - 1][subtree_index_2],
					trees[DEGREE - 1][MAX_RADIUS - 1][subtree_index_2].children_types[nb2],
					DEGREE - 2, MAX_RADIUS - 1);

			subtree_index_3 = trees[DEGREE][MAX_RADIUS][i].cut_off_type;
			subtree_index_3 =
				kill_branch(trees[DEGREE][MAX_RADIUS - 1][subtree_index_3],
				trees[DEGREE - 1][MAX_RADIUS - 1][trees[DEGREE][MAX_RADIUS][i].children_types[nb1]].cut_off_type,
				DEGREE - 1, MAX_RADIUS - 1);
			subtree_index_3 = trees[DEGREE - 1][MAX_RADIUS - 1][subtree_index_3].children_types[nb2];

			//subtree_index_3 = trees[DEGREE][MAX_RADIUS - 1][].children_types[(nb1 <= nb2 ? nb2 + 1 : nb2)];
			//subtree_index_3 = trees[DEGREE - 1][MAX_RADIUS - 1][subtree_index_3].cut_off_type;
			row_index = find_path_tree(subtree_index_1, subtree_index_2, subtree_index_3);

			p = std::find(ia + start_counter, ia + counter, 3 * row_index + 2);
			if(p == ia + counter) {
				write_to_matrix(3 * row_index + 2, i + 1, (-1) * PROB, counter);
			} else {
				// more pointer arithmetic magic
				if(ja[p - ia] != i + 1) { std::cerr << "Error arithmetic magic! Path31\n"; exit(1); }
				ar[p - ia] -= PROB;
			}
			p = std::find(ia + start_counter, ia + counter, 3 * row_index + 3);
			if(p == ia + counter) {
				write_to_matrix(3 * row_index + 3, i + 1, PROB, counter);
			} else {
				if(ja[p - ia] != i + 1) { std::cerr << "Error arithmetic magic! Path32\n"; exit(1); }
				ar[p - ia] += PROB;
			}
		}
	}
	check_epsilon_ar(counter - 1);
	sum_to_one(counter, 3 * p_trees[MAX_RADIUS - 1].size() + 1,
		trees[DEGREE][MAX_RADIUS].size());
	counter--;
	std::cout << "Path polyhedron: Generated " << counter << " coefficients.\n";
	return counter;
}

void print_lp(glp_prob * lp, int row_count, int column_count) {
	int length;
	int ind[column_count + 1];
	double val[column_count + 1];
	std::cout << "\n\nPrinting LP " << glp_get_prob_name(lp) << ":\n\n";
	for(int i = 1; i <= row_count; i++) {
		length = glp_get_mat_row(lp, i, ind, val);
		std::cout << "Row " << i << ": 0 = ";
		for(int j = 1; j <= length; j++) {
			std::cout << val[j] << "*x_" << ind[j];
			if(j < length) std::cout << " + ";
		}
		std::cout << "\n\n";
	}
}

void print_tree(const std::vector<Tree>& trees, int index) {
	const Tree& t = trees[index];
	std::cout << "Tree " << index << " ----\nradius: " << t.radius <<
		"  root_color: " << t.root_color << "  cut_off_type: " <<
		t.cut_off_type << "  children_types: ";
	for(size_t i = 0; i < t.children_types.size(); i++)
		std::cout << t.children_types[i] << " ";
	std::cout << '\n';
}

void print_trees() {
	std::cout << "----------------------------\nRegular trees:\n";
	for(int r = 0; r <= MAX_RADIUS - 1; r++) {
		for(size_t i = 0; i < trees[DEGREE][r].size(); i++)
			print_tree(trees[DEGREE][r], i);
	}
	std::cout << "----------------------------\nIrregular trees:\n";
	for(int r = 0; r <= MAX_RADIUS - 1; r++) {
		for(size_t i = 0; i < trees[DEGREE - 1][r].size(); i++)
			print_tree(trees[DEGREE - 1][r], i);
	}

}

// row_count is number of rows of lp_other
// column_count is number of columns of both matrices
void run_LPs(glp_prob* lp, glp_prob* lp_other, int row_count, int column_count) {
	int counter = 1;
	int length, counter_nonzero;
	int ind[column_count + 1];
	double val[column_count + 1], prob_shortest[column_count + 1];
	int success_shortest = -1;
	int min_nonzero_coefs = 100000;
	double obj;
	// simplex parameter initialization:
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_ERR;

	for(size_t i = 1; i <= row_count; i++) {
		length = glp_get_mat_row(lp_other, i, ind, val);
		std::cout << "LP run " << counter << '\n';
		//std::cout << "Objective function description: ";
		for(int j = 1; j <= length; j++) {
			glp_set_obj_coef(lp, ind[j], val[j]);
			//std::cout << val[j] << "*x_" << ind[j] << (j == length ? ";\n" : " + ");
		}
		for(int b = 0; b < 2; b++) {
			if(b) glp_set_obj_dir(lp, GLP_MIN);
			else glp_set_obj_dir(lp, GLP_MAX);
			glp_simplex(lp, &parm);
			if(glp_get_obj_val(lp) > EPSILON || glp_get_obj_val(lp) < (-EPSILON)) {
				std::cout << "\nSuccess on run " << counter << ".\n";
				std::cout << "Objective function value: " << glp_get_obj_val(lp) << ".\n";
				counter_nonzero = 0;
				for(int j = 1; j <= column_count; j++) {
					obj = glp_get_col_prim(lp, j);
					if(!check_epsilon(obj)) { 
						std::cout << "Graph " << j << ": " << obj << "\n"; 
						counter_nonzero++;
					}
				}
				std::cout << "Nonzero coefficients count: " << counter_nonzero << ".\n";
				//exit(0);
				if(counter_nonzero < min_nonzero_coefs) {
					// I found better solution, i.e., with less nonzero coefficients
					min_nonzero_coefs = counter_nonzero;
					success_shortest = counter;
					for(int j = 1; j <= column_count; j++) {
						prob_shortest[j] = glp_get_col_prim(lp,j);
					}
				}
			}
			counter++;
		}
		for(int j = 1; j <= length; j++) {
			glp_set_obj_coef(lp, ind[j], 0);
		}

	}
	std::cout << "\nMinimum nonzero coefficients count: " << min_nonzero_coefs << ".\n";
	std::cout << "Happened on run " << success_shortest << ".\n";
	if(success_shortest != -1) for(int j = 1; j <= column_count; j++) {
		if(!check_epsilon(prob_shortest[j]))
			std::cout << "Graph " << j << ": " << prob_shortest[j] << '\n';
	}
}

void init_LP(glp_prob* lp1, glp_prob* lp2) {
	glp_set_prob_name(lp1, "A-L polytope");
	glp_set_prob_name(lp2, "Path polytope");

	// I need an additional row for the condition that sum of probs is one
	glp_add_rows(lp1, e_trees[MAX_RADIUS - 1].size() + 1);
	glp_add_rows(lp2, 3 * p_trees[MAX_RADIUS - 1].size() + 1);

	for(size_t i = 1; i <= e_trees[MAX_RADIUS - 1].size(); i++)
		glp_set_row_bnds(lp1, i, GLP_FX, 0.0, 0.0);
	for(size_t i = 1; i <= 3 * p_trees[MAX_RADIUS - 1].size(); i++)
		glp_set_row_bnds(lp2, i, GLP_FX, 0.0, 0.0);

	glp_set_row_name(lp1, e_trees[MAX_RADIUS - 1].size() + 1, "Probability condition");
	glp_set_row_name(lp2, 3 * p_trees[MAX_RADIUS - 1].size() + 1, "Probability condition");
	glp_set_row_bnds(lp1, e_trees[MAX_RADIUS - 1].size() + 1, GLP_FX, 1.0, 1.0);
	glp_set_row_bnds(lp2, 3 * p_trees[MAX_RADIUS - 1].size() + 1, GLP_FX, 1.0, 1.0);

	glp_add_cols(lp1, trees[DEGREE][MAX_RADIUS].size());
	glp_add_cols(lp2, trees[DEGREE][MAX_RADIUS].size());
	for(size_t i = 1; i <= trees[DEGREE][MAX_RADIUS].size(); i++) {
		glp_set_col_name(lp1, i, std::to_string(i).c_str());
		glp_set_col_name(lp2, i, std::to_string(i).c_str());
		glp_set_col_bnds(lp1, i, GLP_DB, 0.0, 1.0);
		glp_set_col_bnds(lp2, i, GLP_DB, 0.0, 1.0);
	}
}

int main(int argc, char* argv[]) {
	generate_graphs(0);
	std::cout << "Generating trees done.\n\n";
	write_to_file();

	glp_prob *lp_AL = glp_create_prob();
	glp_prob *lp_path = glp_create_prob();
	init_LP(lp_AL, lp_path);

	int nonzero_coef_count;
	// generate rerooting polyhedron
	nonzero_coef_count = generate_constraints_AL();
	glp_load_matrix(lp_AL, nonzero_coef_count, ia, ja, ar);
	// clear the matrix
	for(int i = 1; i <= nonzero_coef_count; i++) { ia[i] = 0; ia[i] = 0; ar[i] = 0; }
	// generate A-L polyhedron
	nonzero_coef_count = generate_constraints_path();

	std::cout << "Both polytopes generated!!!\n\n";
	glp_load_matrix(lp_path, nonzero_coef_count, ia, ja, ar);
	print_lp(lp_AL, e_trees[MAX_RADIUS - 1].size() + 1,
		trees[DEGREE][MAX_RADIUS].size());
	print_lp(lp_path, 3 * p_trees[MAX_RADIUS - 1].size() + 1,
		trees[DEGREE][MAX_RADIUS].size());

	//run_LPs(lp_AL, lp_path, 3 * p_trees[MAX_RADIUS - 1].size(), trees[DEGREE][MAX_RADIUS].size());
	run_LPs(lp_path, lp_AL, e_trees[MAX_RADIUS - 1].size(),
		trees[DEGREE][MAX_RADIUS].size());

	glp_delete_prob(lp_AL);
	glp_delete_prob(lp_path);
	return 0;
}
