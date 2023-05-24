#ifndef _KDEFECTIVE_CLIQUE_BB_
#define _KDEFECTIVE_CLIQUE_BB_

#include "Utility.h"
#include "Timer.h"

class kDefectiveClique_BB {
private:
	ui n;
	ui *pstart;
	ui *pstart_R;
	ui *pend;
	ui *edges;
    ui *removed_level;

    ui *degree;
    ui *degree_in_S;

    ui K, UB;
    ui *best_solution;
    ui best_solution_size;

    ui *SR; // union of S and R, where S is at the front
    ui *SR_rid; // reverse ID for SR
    std::queue<ui> Qv;

    ui *neighbors;
    ui *nonneighbors;
    ui *buf;
    ui *buf1;
    ui *buf2;
    char *vis;

public:
    kDefectiveClique_BB() {
    	n = 0;
    	pstart = pstart_R = NULL;
    	pend = NULL;
    	edges = NULL;
        degree = degree_in_S = NULL;
        
        best_solution = NULL;
        K = best_solution_size = UB = 0;

        SR = SR_rid = NULL;
        removed_level = NULL;

        neighbors = nonneighbors = NULL;
        buf = buf1 = buf2 = NULL;
        vis = NULL;
    }

    ~kDefectiveClique_BB() {
    	if(pstart != NULL) {
    		delete[] pstart;
    		pstart = NULL;
    	}
    	if(pstart_R != NULL) {
    		delete[] pstart_R;
    		pstart_R = NULL;
    	}
    	if(pend != NULL) {
    		delete[] pend;
    		pend = NULL;
    	}
    	if(edges != NULL) {
    		delete[] edges;
    		edges = NULL;
    	}
        if(degree != NULL) {
            delete[] degree;
            degree = NULL;
        }
        if(degree_in_S != NULL) {
            delete[] degree_in_S;
            degree_in_S = NULL;
        }
        if(best_solution != NULL) {
            delete[] best_solution;
            best_solution = NULL;
        }
        if(SR != NULL) {
            delete[] SR;
            SR = NULL;
        }
        if(SR_rid != NULL) {
            delete[] SR_rid;
            SR_rid = NULL;
        }
        if(removed_level != NULL) {
        	delete[] removed_level;
        	removed_level = NULL;
        }
        if(neighbors != NULL) {
        	delete[] neighbors;
        	neighbors = NULL;
        }
        if(nonneighbors != NULL) {
        	delete[] nonneighbors;
        	nonneighbors = NULL;
        }
        if(buf != NULL) {
        	delete[] buf;
        	buf = NULL;
        }
        if(buf1 != NULL) {
        	delete[] buf1;
        	buf1 = NULL;
        }
        if(buf2 != NULL) {
        	delete[] buf2;
        	buf2 = NULL;
        }
        if(vis != NULL) {
        	delete[] vis;
        	vis = NULL;
        }
    }

    void load_graph(ui _n, ui *_pstart, ui *_pend, ui *_edges) {
    	n = _n;
        ui m = 0;
        for(ui i = 0;i < n;i ++) m += _pend[i] - _pstart[i];

        assert(pstart == NULL);
        pstart = new ui[n+1];
        edges = new ui[m];

        m = 0;
        for(ui i = 0;i < n;i ++) {
        	pstart[i] = m;
        	for(ui j = _pstart[i];j < _pend[i];j ++) edges[m ++] = _edges[j];
        }
        pstart[n] = m;

        printf("load graph of size n=%u, m=%u (undirected), density=%.5lf\n", n, m/2, double(m)/n/(n-1));
    }

    void kDefectiveClique(ui _K, ui _UB, std::vector<ui> &kDC) {
        K = _K;
        UB = _UB;
        if(K == 0) {
        	printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
        	return ;
        }
        best_solution_size = kDC.size();

        assert(degree == NULL);
		degree = new ui[n+1];
		degree_in_S = new ui[n+1];
		best_solution = new ui[n+1];
		SR = new ui[n+1];
		SR_rid = new ui[n+1];
		removed_level = new ui[n+1];
		pstart_R = new ui[n+1];
		pend = new ui[n+1];
		if(n > K) {
			neighbors = new ui[n+1];
			nonneighbors = new ui[n+1];
			buf = new ui[n+1];
			buf1 = new ui[n+1];
			buf2 = new ui[n+1];
		}
		else {
			neighbors = new ui[K+1];
			nonneighbors = new ui[K+1];
			buf = new ui[K+1];
			buf1 = new ui[K+1];
			buf2 = new ui[K+1];
		}
		vis = new char[n+1];

		memset(vis, 0, sizeof(char)*(n+1));
		memset(degree_in_S, 0, sizeof(ui)*(n+1));
		for(ui i = 0;i < n;i ++) removed_level[i] = n;
		for(ui i = 0;i < n;i ++) SR[i] = SR_rid[i] = i;
		for(ui i = 0;i < n;i ++) {
			pstart_R[i] = pstart[i];
			pend[i] = pstart[i+1];
		}
		while(!Qv.empty()) Qv.pop();

		for(ui i = 0;i < n;i ++) {
			degree[i] = pstart[i+1] - pstart[i];
			// printf("degree %u: %u\n", i, degree[i]);
			if(degree[i] + K < best_solution_size) {
				removed_level[i] = 0;
				Qv.push(i);
			}
		}

		ui R_end = n;
		remove_vertices_and_prune(0, R_end, 0);

		if(R_end != n) printf("Initially pruned %u vertices!\n", n - R_end);

        BB_search(0, R_end, 1);

        if(best_solution_size > kDC.size()) {
            kDC.clear();
            for(int i = 0;i < best_solution_size;i ++) kDC.push_back(best_solution[i]);
        }
    }

    int main(int argc, char *argv[]) {
	    if(argc < 3) {
		    printf("Usage: [1]exe [2]k [3]dir [4 option] lb_of_max_kDefectiveClique_size\n");
		    return 0;
	    }
        K = atoi(argv[1]);
        UB = 100000000;
        if(K == 0) {
        	printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
        	return 0;
        }
        std::vector<ui> kDC;
        if(argc >= 4) {
        	best_solution_size = atoi(argv[3]);
        	printf("initial lb: %u\n", best_solution_size);
        	for(ui i = 0;i < best_solution_size;i ++) kDC.pb(0);
        }
        readGraph_binary(argv[2]);
        printf("Finish reading graph\n");
        Timer t;
        kDefectiveClique(K, UB, kDC);
        printf("Maximum %u-DefectiveClique size: %u, time excluding reading: %s (micro seconds)\n", K, best_solution_size, Utility::integer_to_string(t.elapsed()).c_str());
        return 0;
    }

private:
    void readGraph_binary(char* dir) {
	    FILE *f = Utility::open_file( (std::string(dir) + std::string("/b_degree.bin")).c_str(), "rb");

	    ui tt;
	    fread(&tt, sizeof(ui), 1, f);
	    if(tt != sizeof(ui)) {
		    printf("sizeof int is different: edge.bin(%u), machine(%lu)\n", tt, sizeof(ui));
		    return ;
	    }
	    ui m;
	    fread(&n, sizeof(ui), 1, f);
	    fread(&m, sizeof(ui), 1, f);
	    printf("n = %u, m = %u\n", n, m/2);

	    ui *degree = new ui[n];
	    fread(degree, sizeof(ui), n, f);
	    fclose(f);

	    f = Utility::open_file( (std::string(dir) + std::string("/b_adj.bin")).c_str(), "rb");
	    pstart = new ui[n+1];
	    edges = new ui[m];
	    pstart[0] = 0;
	    for(ui i = 0;i < n;i ++) {
	    	pstart[i+1] = pstart[i] + degree[i];
	    	if(degree[i] > 0) fread(edges+pstart[i], sizeof(ui), degree[i], f);
		}
	    fclose(f);

	    delete[] degree;
    }

    void initialization(ui S_end, ui R_end, ui level) { // only reorganizes neighbors of vertices of R
    	assert(level > 0);
    	while(!Qv.empty()) Qv.pop();
    	ui start_idx = S_end;
    	if(start_idx > 0) -- start_idx;
    	//print_array("SR", SR, 0, R_end, 0);
    	for(ui i = start_idx;i < R_end;i ++) {
    		ui u = SR[i];
    		ui non_neighbors_n = 0, end = pstart_R[u];
    		//if(u == 0) print_neighbors(u, pstart, pend, edges);
    		for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level-1;j ++) {
    			assert(removed_level[edges[j]] == level-1||removed_level[edges[j]] == n);
    			if(removed_level[edges[j]] >= level) {
    				edges[end ++] = edges[j];
    				if(SR_rid[edges[j]] < S_end) {
    					std::swap(edges[pstart_R[u]], edges[end-1]);
    					++ pstart_R[u];
    				}
    			}
    			else nonneighbors[non_neighbors_n ++] = edges[j];
    		}
#ifndef NDEBUG
    		if(degree[u] != end - pstart[u]) {
    			printf("u: %u, pstart[u+1]-pstart[u]: %u, degree[u]: %u, end - pstart[u]: %u\n", u, pstart[u+1]-pstart[u], degree[u], end - pstart[u]);
    		}
    		assert(degree[u] == end-pstart[u]);
#endif
    		for(ui j = 0;j < non_neighbors_n;j ++) edges[end ++] = nonneighbors[j];
    		assert((end < pend[u]&&removed_level[edges[end]] < level-1)||end == pend[u]);
#ifndef NDEBUG
    		for(ui j = end;j < pend[u];j ++) {
    			if(removed_level[edges[j]] >= level) printf("removed_level[edges[j]]: %u, level: %u\n", removed_level[edges[j]], level);
    			assert(removed_level[edges[j]] < level);
    		}
#endif
    	}
    }

    void compute_a_heuristic_solution_and_prune(ui &R_end, ui level) {
    	// the following computes the degeneracy ordering and a heuristic solution
#ifndef NDEBUG
    	for(ui i = 0;i < R_end;i ++) {
    		ui u = SR[i];
    		assert(degree[u] + K >= best_solution_size);
    		assert(pstart[u] == pstart_R[u]);
    		ui end = pstart[u];
    		while(end < pend[u]&&removed_level[edges[end]] >= level) ++ end;
    		for(ui j = end;j < pend[u];j ++) {
    			if(removed_level[edges[j]] >= level) printf("removed_level[edges[j]]: %u, level: %u\n", removed_level[edges[j]], level);
    			assert(removed_level[edges[j]] < level);
    		}
    		if(degree[u] != end - pstart[u]) printf("degree[u]: %u, %u\n", degree[u], end-pstart[u]);
    		assert(degree[u] == end - pstart[u]);
    	}
#endif
		ui *core = neighbors;
		ui *rid = nonneighbors;
		ui *id = buf;
		ui *t_degree = buf1;
		ui total_edges = 0;
		for(ui i = 0;i < R_end;i ++) {
			id[i] = 0;
			t_degree[SR[i]] = degree[SR[i]];
			assert(t_degree[SR[i]] < R_end);
			total_edges += degree[SR[i]];
		}
		for(ui i = 0;i < R_end;i ++) ++ id[t_degree[SR[i]]];
		for(ui i = 1;i < R_end;i ++) id[i] += id[i-1];

		for(ui i = 0;i < R_end;i ++) rid[SR[i]] = -- id[t_degree[SR[i]]];
		for(ui i = 0;i < R_end;i ++) id[rid[SR[i]]] = SR[i];

		ui *degree_start = buf2;
		for(ui i = 0, j = 0;i <= R_end;i ++) {
			while(j < R_end&&t_degree[id[j]] < i) ++ j;
			degree_start[i] = j;
		}

		ui max_core = 0;
		for(ui i = 0;i < R_end;i ++) {
			ui u = id[i];
			assert(degree_start[t_degree[u]] == i);
			if(t_degree[u] > max_core) max_core = t_degree[u];
			core[u] = max_core;

			long long t_n = R_end - i;

			if(t_n*(t_n-1)/2 <= total_edges/2 + K&&R_end - i > best_solution_size) {
				best_solution_size = R_end - i;
				for(ui j = i;j < R_end;j ++) best_solution[j-i] = id[j];
				printf("Degen find a solution of size %u\n", best_solution_size);
			}

			++ degree_start[t_degree[u]];
			if(t_degree[u] == 0) continue;

			degree_start[t_degree[u]-1] = degree_start[t_degree[u]];
			for(ui j = pstart[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) if(rid[edges[j]] > i) {
				ui v = edges[j];
				ui pos1 = degree_start[t_degree[v]], pos2 = rid[v];
				std::swap(id[pos1], id[pos2]);
				rid[id[pos1]] = pos1; rid[id[pos2]] = pos2;
				++ degree_start[t_degree[v]];
				-- t_degree[v];
				total_edges -= 2;
			}
		}

		assert(Qv.empty());
		for(ui i = 0;i < R_end;i ++) if(core[SR[i]] + K < best_solution_size) {
			assert(removed_level[SR[i]] > level);
			removed_level[SR[i]] = level;
			Qv.push(SR[i]);
		}
    }

    ui degeneracy_ordering_and_coloring_adj(ui S_end, ui R_end, ui level, ui *color) {
		ui *rid = buf;
		ui *id = buf1;
		ui *t_degree = color;
		ui max_degree = 0;
		for(ui i = S_end;i < R_end;i ++) {
			ui &d = t_degree[SR[i]] = 0;
			for(ui j = pstart_R[SR[i]];j < pend[SR[i]]&&removed_level[SR[i]] >= level;j ++) {
				assert(SR_rid[SR[i]] >= S_end);
				if(SR_rid[SR[i]] < R_end) ++ d;
			}
			if(d > max_degree) max_degree = d;
		}
		memset(id, 0, sizeof(ui)*(max_degree+1));
		for(ui i = S_end;i < R_end;i ++) ++ id[t_degree[SR[i]]];
		for(ui i = 1;i <= max_degree;i ++) id[i] += id[i-1];

		for(ui i = S_end;i < R_end;i ++) rid[SR[i]] = -- id[t_degree[SR[i]]];
		for(ui i = S_end;i < R_end;i ++) id[rid[SR[i]]] = SR[i];

		ui *degree_start = buf2;
		for(ui i = 0, j = 0;i <= max_degree;i ++) {
			while(j < R_end&&t_degree[id[j]] < i) ++ j;
			degree_start[i] = j;
		}

		ui max_core = 0;
		for(ui i = 0;i < R_end-S_end;i ++) {
			ui u = id[i];
			assert(degree_start[t_degree[u]] == i);
			if(t_degree[u] > max_core) max_core = t_degree[u];

			++ degree_start[t_degree[u]];
			if(t_degree[u] == 0) continue;

			degree_start[t_degree[u]-1] = degree_start[t_degree[u]];
			for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) if(SR_rid[edges[j]] < R_end&&rid[edges[j]] > i) {
				ui v = edges[j];
				ui pos1 = degree_start[t_degree[v]], pos2 = rid[v];
				std::swap(id[pos1], id[pos2]);
				rid[id[pos1]] = pos1; rid[id[pos2]] = pos2;
				++ degree_start[t_degree[v]];
				-- t_degree[v];
			}
		}

    	ui max_color = 0;
    	for(ui i = R_end-S_end;i > 0;i --) {
    		ui u = id[i-1];
    		for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) {
    			if(SR_rid[edges[j]] < R_end&&rid[edges[j]] >= i) vis[color[edges[j]]] = 1;
    		}
    		for(ui j = 0;;j ++) if(!vis[j]) {
    			color[u] = j;
    			if(j > max_color) max_color = j;
    			break;
    		}
    		for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) {
    			if(SR_rid[edges[j]] < R_end&&rid[edges[j]] >= i) vis[color[edges[j]]] = 0;
    		}
    	}

    	return max_color + 1;
    }

    bool is_kDefectiveClique(ui R_end) {
    	ui total_edges = 0;
    	for(ui i = 0;i < R_end;i ++) total_edges += degree[SR[i]];
    	long long all_edges = R_end;
    	return all_edges*(R_end-1)/2 <= total_edges/2 + K;
    }

    void store_a_kDefectiveClique(ui S_end) {
    	assert(S_end > best_solution_size);
		best_solution_size = S_end;
		for(ui i = 0;i < best_solution_size;i ++) best_solution[i] = SR[i];
#ifndef NDEBUG
		printf("Find a kDefectiveClique of size: %u\n", best_solution_size);
#endif
    }

    void BB_search(ui S_end, ui R_end, ui level) {
#ifndef NDEBUG
    	assert(compute_missing_edges_in_S(S_end) <= K);
    	for(ui i = 0;i < R_end;i ++) assert(degree[SR[i]] + K >= best_solution_size);
    	//for(ui i = 0;i < n;i ++) assert(!vis[i]);
#endif

        if(best_solution_size >= UB||R_end <= best_solution_size) return ;

        fflush(stdout);

    	if(S_end > best_solution_size) store_a_kDefectiveClique(S_end);
        if(R_end > best_solution_size&&is_kDefectiveClique(R_end)) store_a_kDefectiveClique(R_end);
        if(R_end <= best_solution_size+1) return ;

        // printf("level %u\n", level);

        initialization(S_end, R_end, level);

#ifndef NDEBUG
        for(ui i = S_end;i < R_end;i ++) {
			ui u = SR[i], cnt = 0;
			for(ui j = pstart[u];j < pstart_R[u];j ++) assert(SR_rid[edges[j]] < S_end);
			assert(degree_in_S[u] == pstart_R[u] - pstart[u]);
			for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) if(SR_rid[edges[j]] < R_end) {
				assert(SR_rid[edges[j]] >= S_end);
				++ cnt;
			}
			assert(degree[u] == cnt + degree_in_S[u]);
		}
#endif

        ui old_R_end = R_end;

        assert(Qv.empty());
        bool terminate = collect_removable_vertices_based_on_degree_in_S(S_end, R_end, level);
        if(terminate||remove_vertices_and_prune(S_end, R_end, level)) {
        	restore_R(S_end, R_end, old_R_end, level);
        	return ;
        }

        //printf("here\n");

        if(S_end == 0) {
        	assert(Qv.empty());
        	compute_a_heuristic_solution_and_prune(R_end, level);
        	if(level == 1) printf("First level kdefectiveclique size: %lu\n", best_solution_size);
        	if(remove_vertices_and_prune(S_end, R_end, level)) {
        		restore_R(S_end, R_end, old_R_end, level);
        		return ;
        	}
        }

        if(S_end > 0) {
        	assert(Qv.empty());
        	collect_removable_vertices_based_on_vertex_pair(S_end, R_end, level); // remove vertices from R
			if(remove_vertices_and_prune(S_end, R_end, level)) {
        		restore_R(S_end, R_end, old_R_end, level);
				return ;
			}
        }

        //printf("here1\n");

		if(R_end > best_solution_size&&is_kDefectiveClique(R_end)) store_a_kDefectiveClique(R_end);
		if(R_end <= best_solution_size+1) {
        	restore_R(S_end, R_end, old_R_end, level);
			return ;
		}

		if(degree_based_prune(S_end, R_end, level)||coloring_based_prune(S_end, R_end, level)) {
			restore_R(S_end, R_end, old_R_end, level);
			return ;
		}

		//printf("here2\n");

        ui u = n; // u is the branching vertex
        bool must_include = false;
        for(ui i = S_end;i < R_end;i ++) if(degree[SR[i]] + 2 >= R_end) {
       		must_include = true;
       		u = SR[i];
       		break;
       	}
		if(u == n) u = choose_branch_vertex_based_on_degree(S_end, R_end, level);
		assert(u != n&&SR_rid[u] >= S_end);
		assert(SR_rid[u] < R_end);
		assert(SR[SR_rid[u]] == u);
        assert(degree[u] + K >= best_solution_size);

#ifndef NDEBUG
        for(ui i = S_end;i < R_end;i ++) {
			ui u = SR[i], cnt = 0;
			for(ui j = pstart[u];j < pstart_R[u];j ++) assert(SR_rid[edges[j]] < S_end);
			assert(degree_in_S[u] == pstart_R[u] - pstart[u]);
			for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) if(SR_rid[edges[j]] < R_end) {
				assert(SR_rid[edges[j]] >= S_end);
				++ cnt;
			}
			assert(degree[u] == cnt + degree_in_S[u]);
		}
    	//for(ui i = 0;i < n;i ++) assert(!vis[i]);
#endif

        //printf("here3\n");

        // the first branch includes u into S
        ui pre_best_solution_size = best_solution_size, t_old_R_end = R_end;
        move_u_from_R_to_S(u, S_end, R_end, level);
        if(!remove_vertices_and_prune(S_end, R_end, level)) BB_search(S_end, R_end, level+1);
        restore_R(S_end, R_end, t_old_R_end, level);

        if(must_include) {
        	move_u_from_S_to_R(S_end, R_end, level);
        	restore_R(S_end, R_end, old_R_end, level);
        	return ;
        }

        //printf("back to level: %u\n", level);
#ifndef NDEBUG
        for(ui i = S_end;i < R_end;i ++) {
			ui u = SR[i], cnt = 0;
			for(ui j = pstart[u];j < pstart_R[u];j ++) assert(SR_rid[edges[j]] < S_end);
			//assert(degree_in_S[u] == pstart_R[u] - pstart[u]);
			for(ui j = pstart_R[u];j < pend[u]&&removed_level[edges[j]] >= level;j ++) if(SR_rid[edges[j]] < R_end) {
				if(SR_rid[edges[j]] >= S_end) ++ cnt;
			}
			assert(degree[u] == cnt + degree_in_S[u]);
		}
    	//for(ui i = 0;i < n;i ++) assert(!vis[i]);
#endif

        // the second branch excludes u from S
        assert(Qv.empty());
        bool pruned = remove_u_from_S_with_prune(S_end, R_end, level);
        if(!pruned&&best_solution_size > pre_best_solution_size) pruned = collect_removable_vertices(S_end, R_end, level);
        if(!pruned) {
            if(!remove_vertices_and_prune(S_end, R_end, level)) BB_search(S_end, R_end, level+1);
        }
        restore_R(S_end, R_end, old_R_end, level);
    }

    bool collect_removable_vertices_based_on_degree_in_S(ui S_end, ui R_end, ui level) {
    	ui *cnt = buf;
    	memset(cnt, 0, sizeof(ui)*(S_end+1));
    	for(ui i = S_end;i < R_end;i ++) ++ cnt[S_end - degree_in_S[SR[i]]];
    	assert(R_end > best_solution_size);

    	ui missing_edges = compute_missing_edges_in_S(S_end);
    	ui remaining_vertices_n = best_solution_size - S_end, idx = S_end;
    	assert(missing_edges <= K&&S_end <= best_solution_size);
    	assert(Qv.empty());
    	for(ui i = 0;i <= S_end;i ++) {
    		if(cnt[i] < remaining_vertices_n) {
    			remaining_vertices_n -= cnt[i];
    			missing_edges += i*cnt[i];
    		}
    		else {
    			missing_edges += i*remaining_vertices_n;
    			idx = i;

    			ui next_value = i;
    			if(cnt[i] == remaining_vertices_n) {
    				for(i ++;i <= S_end;i ++) if(cnt[i]) {
    					next_value = i;
    					break;
    				}
    			}
    			if(missing_edges + next_value > K) return true;

    			break;
    		}
    	}
    	for(ui i = S_end;i < R_end;i ++) if(S_end-degree_in_S[SR[i]] > idx&&missing_edges+S_end-degree_in_S[SR[i]] > K) {
    		removed_level[SR[i]] = level;
    		Qv.push(SR[i]);
    	}

    	return false;
    }

    void collect_removable_vertices_based_on_vertex_pair(ui S_end, ui R_end, ui level) {
       	assert(S_end >= 1&&Qv.empty());
       	ui u = SR[S_end-1], neighbors_n = 0, non_neighbors_n = 0;
       	get_neighbors_and_non_neighbors_in_R(u, S_end, R_end, level, neighbors_n, non_neighbors_n);
       	for(ui i = 0;i < neighbors_n;i ++) vis[neighbors[i]] = 1;
       	for(ui i = 0;i < non_neighbors_n;i ++) vis[nonneighbors[i]] = 2;

       	ui missing_edges_in_S = compute_missing_edges_in_S(S_end);
       	assert(missing_edges_in_S <= K);

       	for(ui i = 0;i < neighbors_n;i ++) {
       		ui v = neighbors[i];
       		ui common_neighbors = 0, exclusive_non_neighbor_u = 0;
       		for(ui j = pstart_R[v];j < pend[v]&&removed_level[edges[j]] >= level;j ++) {
       			if(vis[edges[j]] == 1) ++ common_neighbors;
       			else if(vis[edges[j]] == 2) ++ exclusive_non_neighbor_u;
       		}
       		ui exclusive_non_neighbor_v = neighbors_n - 1 - common_neighbors;
       		ui exclusive_non_neighbors = exclusive_non_neighbor_u + exclusive_non_neighbor_v;
       		ui common_non_neighbors = non_neighbors_n - exclusive_non_neighbor_u;

       		ui UB = S_end + 1 + common_neighbors + mmin(K-missing_edges_in_S, exclusive_non_neighbors);
       		if(exclusive_non_neighbors < K-missing_edges_in_S) {
       			ui tmp = (K-missing_edges_in_S-exclusive_non_neighbors)/2;
       			UB += mmin(tmp, common_non_neighbors);
       		}
       		if(UB <= best_solution_size) {
       			removed_level[v] = level;
       			vis[v] = 0;
       			Qv.push(v);
       			neighbors[i] = neighbors[-- neighbors_n];
       			-- i;
       		}
       	}

       	++ missing_edges_in_S;
       	assert(missing_edges_in_S <= K||non_neighbors_n == 0);
       	for(ui i = 0;i < non_neighbors_n;i ++) {
   			ui v = nonneighbors[i];
   			ui common_neighbors = 0, exclusive_non_neighbor_u = 0;
   			for(ui j = pstart[v];j < pend[v]&&removed_level[edges[j]] >= level;j ++) {
				if(vis[edges[j]] == 1) ++ common_neighbors;
				else if(vis[edges[j]] == 2) ++ exclusive_non_neighbor_u;
			}
			ui exclusive_non_neighbor_v = neighbors_n - common_neighbors;
       		ui exclusive_non_neighbors = exclusive_non_neighbor_u + exclusive_non_neighbor_v;
			ui common_non_neighbors = non_neighbors_n - 1 - exclusive_non_neighbor_u;

       		ui UB = S_end + 1 + common_neighbors + mmin(K-missing_edges_in_S, exclusive_non_neighbors);
       		if(exclusive_non_neighbors < K-missing_edges_in_S) {
       			ui tmp = (K-missing_edges_in_S-exclusive_non_neighbors)/2;
       			UB += mmin(tmp, common_non_neighbors);
       		}
   			if(UB <= best_solution_size) {
   				removed_level[v] = level;
   				vis[v] = 0;
   				Qv.push(v);
   				nonneighbors[i] = nonneighbors[-- non_neighbors_n];
   				-- i;
   			}
   		}
       	for(ui i = 0;i < neighbors_n;i ++) vis[neighbors[i]] = 0;
       	for(ui i = 0;i < non_neighbors_n;i ++) vis[nonneighbors[i]] = 0;
    }

    bool collect_removable_vertices(ui S_end, ui R_end, ui level) {
		for(ui i = 0;i < S_end;i ++) if(degree[SR[i]] + K < best_solution_size) return true;

		for(ui i = S_end;i < R_end;i ++) if(removed_level[SR[i]] > level){
			ui v = SR[i];
			if(degree[v] + K < best_solution_size) {
				removed_level[v] = level;
				Qv.push(v);
			}
		}

		return false;
	}

    bool degree_based_prune(ui S_end, ui R_end, ui level) {
    	ui *cnt = buf;
    	memset(cnt, 0, sizeof(ui)*(S_end+1));
    	for(ui i = S_end;i < R_end;i ++) ++ cnt[S_end - degree_in_S[SR[i]]];

    	ui missing_edges = compute_missing_edges_in_S(S_end);
    	ui t_UB = S_end;
    	for(ui i = 0;i <= S_end;i ++) for(ui j = 0;j < cnt[i];j ++) {
    		if(missing_edges + i > K) break;

    		missing_edges += i;
    		++ t_UB;
    	}

    	return t_UB <= best_solution_size;
    }

    bool coloring_based_prune(ui S_end, ui R_end, ui level) {
#ifndef NDEBUG
    	//for(ui i = 0;i < n;i ++) assert(!vis[i]);
#endif
    	ui *color = neighbors;
    	ui color_n = degeneracy_ordering_and_coloring_adj(S_end, R_end, level, color);
    	// printf("color_n: %u\n", color_n);


    	ui missing_edges_n = compute_missing_edges_in_S(S_end);
    	assert(missing_edges_n <= K);
    	if(S_end + color_n + K - missing_edges_n <= best_solution_size) return true;

#ifndef NDEBUG
    	// for(ui i = 0;i < n;i ++) assert(!vis[i]);
    	for(ui i = S_end;i < R_end;i ++) vis[color[SR[i]]] = 1;
    	for(ui i = 0;i < color_n;i ++) assert(vis[i]);
    	for(ui i = S_end;i < R_end;i ++) vis[color[SR[i]]] = 0;
    	//for(ui i = S_end;i < R_end;i ++) if(color[SR[i]] == 16) printf("color of %u is 16\n", SR[i]);
#endif

    	ui *head = nonneighbors;
    	ui *next = buf;
    	for(ui i = 0;i <= S_end;i ++) head[i] = n;
    	for(ui i = S_end;i < R_end;i ++) {
    		assert(degree_in_S[SR[i]] <= S_end);
    		ui u = SR[i], non_neighbors_n = S_end - degree_in_S[SR[i]];
    		//if(color[u] == 16) printf("color[%u] is 16\n", u);
    		next[u] = head[non_neighbors_n];
    		head[non_neighbors_n] = u;
    	}

    	ui *color_head = buf1;
    	ui *color_next = buf2;
    	for(ui i = 0;i < color_n;i ++) color_head[i] = n;
    	for(ui i = 0;i <= S_end;i ++) for(ui u = head[S_end-i];u != n;u = next[u]) {
    		ui c = color[u];
    		assert(c < color_n);
    		color_next[u] = color_head[c];
    		color_head[c] = u;
    		//printf("color_head[%u] = %u\n", c, u);
    	}

    	for(ui i = 0;i <= K;i ++) head[i] = n;
    	for(ui i = 0;i < color_n;i ++) {
    		ui u = color_head[i];
    		//printf("color_head[%u]: %u\n", i, u);
    		assert(u != n&&color[u] == i);
    		color_head[i] = color_next[u];

    		ui non_neighbors_n = S_end - degree_in_S[u];
    		assert(non_neighbors_n <= K);
    		next[u] = head[non_neighbors_n];
    		head[non_neighbors_n] = u;
    	}

    	ui min_key = 0, UB = S_end;
    	while(true) {
    		while(min_key <= K&&head[min_key] == n) ++ min_key;
    		if(min_key > K||missing_edges_n + min_key > K) break;

    		missing_edges_n += min_key;
    		++ UB;
    		ui u = head[min_key];
    		assert(u != n);
    		head[min_key] = next[u];
    		ui new_u = color_head[color[u]];
    		if(new_u == n) continue;

    		color_head[color[u]] = color_next[new_u];
    		ui non_neighbors_n = S_end - degree_in_S[u];
    		ui new_key = S_end - degree_in_S[new_u] + 1 + min_key - non_neighbors_n;
    		next[new_u] = head[new_key];
    		head[new_key] = new_u;
    	}

    	return UB <= best_solution_size;
    }

    void get_neighbors_and_non_neighbors_in_R(ui u, ui S_end, ui R_end, ui level, ui &neighbors_n, ui &non_neighbors_n) {
    	neighbors_n = non_neighbors_n = 0;
    	for(ui i = pstart_R[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++) if(SR_rid[edges[i]] < R_end) vis[edges[i]] = 1;
    	for(ui i = S_end;i < R_end;i ++) if(SR[i] != u) {
    		if(vis[SR[i]]) neighbors[neighbors_n++] = SR[i];
    		else nonneighbors[non_neighbors_n++] = SR[i];
    	}
    	for(ui i = pstart_R[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++) vis[edges[i]] = 0;
    }

    void move_u_from_R_to_S(ui u, ui &S_end, ui R_end, ui level) {
    	assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end&&SR[SR_rid[u]] == u);
    	swap_pos(S_end, SR_rid[u]);
    	++ S_end;

    	for(ui i = pstart[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++) if(SR_rid[edges[i]] < R_end) {
    		++ degree_in_S[edges[i]];
    	}

        ui missing_edges = compute_missing_edges_in_S(S_end);
        assert(missing_edges <= K);

        assert(Qv.empty());
        for(ui i = S_end;i < R_end;i ++) if(S_end - degree_in_S[SR[i]] + missing_edges > K) {
        	removed_level[SR[i]] = level;
        	Qv.push(SR[i]);
        }
    }

    void move_u_from_S_to_R(ui &S_end, ui R_end, ui level) {
		assert(S_end > 0);
		ui u = SR[-- S_end];
		for(ui i = pstart[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++) if(SR_rid[edges[i]] < R_end) {
			ui v = edges[i];
			-- degree_in_S[v];
			if(SR_rid[v] >= S_end&&pstart_R[v] > pstart[v]&&edges[pstart_R[v]-1] == u) -- pstart_R[v];
		}
	}

    bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level) {
		assert(S_end > 0);
		ui u = SR[S_end-1];
		-- S_end; -- R_end;
		swap_pos(S_end, R_end);
		removed_level[u] = level;

		bool ret = false;
		for(ui i = pstart[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++) if(SR_rid[edges[i]] < R_end) {
			ui v = edges[i];
			-- degree_in_S[v];
			if(SR_rid[v] >= S_end&&pstart_R[v] > pstart[v]&&edges[pstart_R[v]-1] == u) -- pstart_R[v];

			-- degree[v];
			if(degree[v] + K < best_solution_size) {
				if(SR_rid[v] < S_end) ret = true;
				else {
					assert(removed_level[v] > level);
					removed_level[v] = level;
					Qv.push(v);
				}
			}
		}
		return ret;
	}

    bool remove_vertices_and_prune(ui S_end, ui &R_end, ui level) {
		while(!Qv.empty()) {
			ui u = Qv.front(); Qv.pop(); // remove u
			assert(SR[SR_rid[u]] == u);
			assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end);
			-- R_end;
			swap_pos(SR_rid[u], R_end);

			// printf("remove %u with degree %u\n", u, degree[u]);

			bool terminate = false;
			for(ui i = pstart[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++) if(SR_rid[edges[i]] < R_end) {
				ui w = edges[i];
				assert(degree[w] > 0);
				-- degree[w];
				if(degree[w] + K < best_solution_size) {
					if(SR_rid[w] < S_end) terminate = true; // UB
					else if(removed_level[w] > level) { // RR
						removed_level[w] = level;
						Qv.push(w);
					}
				}
			}
			if(terminate) return true;
		}

        return false;
    }

    void restore_R(ui S_end, ui &R_end, ui old_R_end, ui level) {
        while(!Qv.empty()) {
            ui u = Qv.front(); Qv.pop();
            assert(removed_level[u] == level&&SR_rid[u] >= S_end&&SR_rid[u] < R_end);
            removed_level[u] = n;
        }
        while(R_end < old_R_end) { // insert u back into R
            ui u = SR[R_end];
            assert(removed_level[u] == level&&SR_rid[u] == R_end);
            removed_level[u] = n;

            for(ui i = pstart[u];i < pend[u]&&removed_level[edges[i]] >= level;i ++) {
            	assert(SR[SR_rid[edges[i]]] == edges[i]);
            	if(SR_rid[edges[i]] < R_end) ++ degree[edges[i]];
            }

            ++ R_end;
            // printf("restore %u with degree %u\n", u, degree[u]);
        }
    }
    
    ui compute_missing_edges_in_S(ui S_end) {
    	ui res = 0;
    	for(ui i = 0;i < S_end;i ++) {
    		assert(degree_in_S[SR[i]] < S_end);
    		res += S_end - 1 - degree_in_S[SR[i]];
    	}
    	return res/2;
    }

    void swap_pos(ui i, ui j) {
        std::swap(SR[i], SR[j]);
        SR_rid[SR[i]] = i;
        SR_rid[SR[j]] = j;
    }

    ui choose_branch_vertex_based_on_degree(ui S_end, ui R_end, ui level) {
    	//printf("start choosing branching vertices\n");
    	ui u = SR[S_end];
		for(ui i = S_end+1;i < R_end;i ++) {
			ui v = SR[i];
			if(degree_in_S[u] == S_end) {
				if(degree_in_S[v] != S_end||degree[v] > degree[u]) u = v;
			}
			else if(degree_in_S[v] != S_end&&(degree[v] > degree[u]||(degree[v] == degree[u]&&degree_in_S[v] < degree_in_S[u]))) u = v;	
		}

#ifndef NDEBUG
		ui missing_edges_n = compute_missing_edges_in_S(S_end);
		assert(u != n&&missing_edges_n + S_end - degree_in_S[u] <= K);
		assert(degree[u] + 2 < R_end);
#endif

		return u;
    }

    void print_neighbors(ui u, const ui *pstart, const ui *pend, const ui *edges) {
    	std::vector<ui> neighbors;
    	for(ui i = pstart[u];i < pend[u];i ++) neighbors.push_back(edges[i]);
    	//sort(neighbors.begin(), neighbors.end());
    	printf("neighbors of %u:", u);
    	for(ui i = 0;i < neighbors.size();i ++) printf(" %u(%u)", neighbors[i], removed_level[neighbors[i]]);
    	printf("\n");
    }

    void print_array(const char *str, const ui *array, ui idx_start, ui idx_end, ui l) {
    	for(ui i = 0;i < l;i ++) printf(" ");
    	printf("%s:", str);
    	for(ui i = idx_start;i < idx_end;i ++) printf(" %u", array[i]);
    	printf("\n");
    }
};

#endif
