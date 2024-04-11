#include <bits/stdc++.h>
using namespace std;
const int tissue_count = 10;
map<long long, int> vertex_to_number; // mapping of vertex to number
map<int, long long> number_to_vertex; // mapping of vertex to number
vector<vector<pair<int, long long>>> graph_with_edges;
char types[tissue_count][8] = {"AD2", "AO", "BL1", "LG", "LI11", "OV2", "PA", "PO3", "SB", "SX"};
const int bitsetsize = 1450000; // max number of vertices in a chromosome
bitset<bitsetsize> interesting_vertices[1 << tissue_count];
vector<tuple<int, long long, int>> result;
int main(int argc, char *argv[])
{

	if (argc < 8)
	{
		printf("Not enough arguments.\n");
		return 0;
	}
	if (argc > 8)
	{
		printf("Too many command line arguments.\n");
		return 0;
	}
	// nmax nmin a b

	char *inputfilename = argv[1];
	char *outputfilename = argv[2];
	char *csv_outputfilename = argv[3];
	int nmax, nmin;
	double a, b;

	nmax = atoi(argv[4]);
	nmin = atoi(argv[5]);
	a = atof(argv[6]);
	b = atof(argv[7]);

	FILE *inputfile;
	FILE *outputfile;
	FILE *csv_outputfile;

	inputfile = fopen(inputfilename, "r");
	outputfile = fopen(outputfilename, "w");

	char buff[255];
	fscanf(inputfile, "%s", buff);

	int chr;
	long long bait, oe;
	set<int> vertices;

	while (fscanf(inputfile, "%d,%lld,%lld,%s", &chr, &bait, &oe, buff) != EOF)
	{
		vertices.insert(bait);
		vertices.insert(oe);
	}

	if (vertices.size() > bitsetsize)
	{
		printf("Vertices don't fit in current bitset size\n");
		return 0;
	}

	int j = 0;
	for (auto i : vertices)
	{
		vertex_to_number[i] = j;
		number_to_vertex[j] = i;
		j++;
		vector<pair<int, long long>> tmp;
		graph_with_edges.push_back(tmp);
	}

	fclose(inputfile);
	inputfile = fopen(inputfilename, "r");
	fscanf(inputfile, "%s", buff);

	while (fscanf(inputfile, "%d,%lld,%lld", &chr, &bait, &oe) != EOF)
	{
		long long mask = 0;
		int tmp;
		for (int i = 0; i < tissue_count; i++)
		{
			fscanf(inputfile, ",%d", &tmp);
			mask *= 2;
			mask += tmp;
		}

		graph_with_edges[vertex_to_number[bait]].push_back(make_pair(vertex_to_number[oe], mask));
		graph_with_edges[vertex_to_number[oe]].push_back(make_pair(vertex_to_number[bait], mask));
	}

	printf("Graph has been read\n");

	csv_outputfile = fopen(csv_outputfilename, "w");
	fprintf(csv_outputfile, "V,E,F,countbasetissues,basetissues,countatissues,atissues,countbtissues,btissues,from,to\n");
	for (long long mask = 1; mask < (1ll << tissue_count); mask++)
	{

		for (int i = 0; i < tissue_count; i++)
		{
			if ((mask & (1 << i)) == (1 << i))
			{
				interesting_vertices[mask] |= interesting_vertices[mask ^ (1 << i)];
			}
		}

		vector<vector<int>> components;
		vector<bool> visited;
		int visited_count = 0;
		for (auto j : vertices)
		{
			visited.push_back(0);
		}

		for (int j = 0; j < vertices.size(); j++)
		{
			if (visited[j] == true || interesting_vertices[mask][j] == true)
			{
				continue;
			}

			visited[j] = true;
			vector<int> vertex_vector;
			vertex_vector.push_back(j);
			vector<int> single_component;
			single_component.push_back(j);
			while (vertex_vector.size() > 0)
			{
				int tmp_vertex = vertex_vector[vertex_vector.size() - 1];
				vertex_vector.pop_back();

				for (auto k : graph_with_edges[tmp_vertex])
				{
					int n_vertice = get<0>(k);

					if (!visited[n_vertice] && (get<1>(k) & mask) == mask)
					{
						vertex_vector.push_back(n_vertice);
						single_component.push_back(n_vertice);
						visited[n_vertice] = true;
					}
				}
			}

			components.push_back(single_component);
		}

		for (auto component : components)
		{
			if (component.size() >= nmin && component.size() <= nmax)
			{

				int component_edge_count = 0;
				for (auto v : component)
				{
					for (auto e : graph_with_edges[v])
					{
						if ((get<1>(e) & mask) == mask)
							component_edge_count++;
					}
				}

				int less_count = 0, more_count = 0;
				int base_count = 0;
				long long full_mask = mask;
				bool less_types[tissue_count], more_types[tissue_count];
				int one_edge_end, another_edge_end;

				/// Get aTissues and bTissues for one baseTissue set
				for (int j = 0; j < tissue_count; j++)
				{
					long long tmpmask = mask | (1ll << j);
					if (tmpmask == mask)
					{
						base_count++;
						more_types[j] = false;
						less_types[j] = false;
						continue;
					}
					int edge_count = 0;
					for (auto v : component)
					{
						for (auto e : graph_with_edges[v])
						{
							if ((get<1>(e) & tmpmask) == tmpmask)
								edge_count++;
						}
					}

					if (edge_count == component_edge_count)
					{
						base_count++;
						full_mask |= (1ll << j);
					}
					else if (edge_count >= component_edge_count * a)
					{
						more_count++;
						more_types[j] = true;
						less_types[j] = false;
					}
					else if (edge_count <= component_edge_count * b)
					{
						less_count++;
						more_types[j] = false;
						less_types[j] = true;
					}
					else
					{
						more_types[j] = false;
						less_types[j] = false;
					}
				}

				if (less_count * (more_count + base_count) > 0)
				{
					for (int i = 0; i < component.size(); i++)
					{
						interesting_vertices[mask][component[i]] = true;
					}

					int tmp = less_count * (more_count + base_count);
					result.push_back(tie(tmp, full_mask, component[0]));
				}
			}
		}
		if (mask % 1000 == 0)
		{
			printf("%.2lf%% processed, %d components found.\n", (double)mask / (double)(1 << tissue_count) * 100, result.size());
		}
	}

	sort(result.begin(), result.end());

	for (int i = result.size() - 1; i >= 0; i--)
	{
		auto r = result[i];
		// Here mask is the full mask (all tissues with 100%)
		long long mask = get<1>(r);
		int vertex = get<2>(r);

		vector<bool> visited;

		for (auto j : vertices)
		{
			visited.push_back(0);
		}

		visited[vertex] = true;
		vector<int> vertex_vector;
		vertex_vector.push_back(vertex);
		vector<int> component;
		component.push_back(vertex);
		while (vertex_vector.size() > 0)
		{
			int tmp_vertex = vertex_vector[vertex_vector.size() - 1];
			vertex_vector.pop_back();

			for (auto k : graph_with_edges[tmp_vertex])
			{
				int n_vertice = get<0>(k);
				if (!visited[n_vertice] && ((get<1>(k) & mask) == mask))
				{
					vertex_vector.push_back(n_vertice);
					component.push_back(n_vertice);
					visited[n_vertice] = true;
				}
			}
		}

		int component_edge_count = 0;
		for (auto v : component)
		{
			for (auto e : graph_with_edges[v])
			{
				if ((get<1>(e) & mask) == mask)
					component_edge_count++;
			}
		}

		int less_count = 0, more_count = 0;
		int base_count = 0;
		bool less_types[tissue_count], more_types[tissue_count];
		int one_edge_end, another_edge_end;

		for (int j = 0; j < tissue_count; j++)
		{
			long long tmpmask = mask | (1ll << j);
			if (tmpmask == mask)
			{
				base_count++;
				more_types[j] = false;
				less_types[j] = false;
			}
			int edge_count = 0;
			for (auto v : component)
			{
				for (auto e : graph_with_edges[v])
				{
					if ((get<1>(e) & tmpmask) == tmpmask)
					{
						edge_count++;
						one_edge_end = v;
						another_edge_end = get<0>(e);
					}
				}
			}
			if (edge_count >= component_edge_count * a && edge_count < component_edge_count)
			{
				more_count++;
				more_types[j] = true;
				less_types[j] = false;
			}
			else if (edge_count <= component_edge_count * b)
			{
				less_count++;
				more_types[j] = false;
				less_types[j] = true;
			}
			else
			{
				more_types[j] = false;
				less_types[j] = false;
			}
		}

		int tmp = (base_count + more_count) * less_count;
		double sqrt_result = sqrt(tmp);

		fprintf(outputfile, "1 %d %d %.2lf\n", component.size(), component_edge_count / 2, sqrt_result); // list of vertices

		fprintf(outputfile, "2 ");
		int basecount = 0;
		long long atissues = 0, btissues = 0;
		long long counter = 1;
		for (int j = 0; j < tissue_count; j++)
		{
			if ((mask & (1 << j)) == (1 << j))
			{
				fprintf(outputfile, "%s ", types[tissue_count - 1 - j]);
				basecount++;
			}
		}
		counter = 1;
		fprintf(outputfile, "\n3 ");
		for (int j = 0; j < tissue_count; j++)
		{
			atissues *= 2;
			if (more_types[j] == true)
			{
				fprintf(outputfile, "%s ", types[tissue_count - 1 - j]);
				atissues += 1;
			}
		}

		fprintf(outputfile, "\n4 ");
		for (int j = 0; j < tissue_count; j++)
		{
			btissues *= 2;
			if (less_types[j] == true)
			{
				fprintf(outputfile, "%s ", types[tissue_count - 1 - j]);
				btissues += 1;
			}
		}
		fprintf(outputfile, "\n");

		fprintf(outputfile, "5 %lld %lld\n", number_to_vertex[one_edge_end], number_to_vertex[another_edge_end]);

		long long amask = 0, bmask = 0;
		long long MASK = 0;
		counter = 1;
		for (int j = 0; j < tissue_count; j++, counter *= 2)
		{
			MASK *= 2;

			if ((mask / counter) % 2 == 1)
			{
				MASK += 1;
			}
		}
		amask = atissues;
		bmask = btissues;

		fprintf(csv_outputfile, "%d,%d,%.2lf,%d,%lld,", component.size(), component_edge_count / 2, sqrt_result, basecount, MASK);
		fprintf(csv_outputfile, "%d,%lld,%d,%lld,%lld,%lld\n", more_count, amask, less_count, bmask, number_to_vertex[one_edge_end], number_to_vertex[another_edge_end]);
	}

	printf("Done\n");
	return 0;
}
