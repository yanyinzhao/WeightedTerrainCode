#include <iostream>
#include <fstream>
#include <string>

#include "geodesic_algorithm_subdivision.h"
#include "geodesic_algorithm_subdivision_log.h"
#include "divide_and_conquer.h"

int main(int argc, char **argv)
{
	if (argc < 6)
	{
		std::cout << "run: ./main [terrain_data] [epsilon_SP] [epsilon_SL] [calculate_exact_path] [calculate_FSP]" << std::endl;
		return 0;
	}

	std::string input_dataset = argv[1];
	double Steiner_point_epsilon = std::stod(argv[2]);
	double snell_law_epsilon = std::stod(argv[3]);
	int calculate_exact_path = std::stoi(argv[4]); // 1 for calculate exact path, 0 for not calculate exact path
	int calculate_FSP = std::stoi(argv[5]);		   // 1 for calculate FSP, 0 for not calculate FSP

	std::vector<double> points;
	std::vector<unsigned> faces;

	std::string input_folder = "../input/";
	std::string input_file = input_folder + input_dataset;

	int source_index;
	int destination_index;

	if (input_dataset.compare("BH.off") == 0)
	{
		source_index = 3139;
		destination_index = 2705;
	}
	else if (input_dataset.compare("EP.off") == 0)
	{
		source_index = 5859;
		destination_index = 6356;
	}
	else if (input_dataset.compare("BH_small.off") == 0)
	{
		source_index = 1007;
		destination_index = 768;
	}
	else if (input_dataset.compare("EP_small.off") == 0)
	{
		source_index = 767;
		destination_index = 1071;
	}

	else if (input_dataset.compare("EP_200192.off") == 0)
	{
		source_index = 14527;
		destination_index = 72453;
	}
	else if (input_dataset.compare("EP_401170.off") == 0)
	{
		source_index = 28645;
		destination_index = 145290;
	}
	else if (input_dataset.compare("EP_600996.off") == 0)
	{
		source_index = 43052;
		destination_index = 217308;
	}
	else if (input_dataset.compare("EP_800768.off") == 0)
	{
		source_index = 57340;
		destination_index = 289612;
	}
	else if (input_dataset.compare("EP_1000768.off") == 0)
	{
		source_index = 71991;
		destination_index = 361379;
	}
	else if (input_dataset.compare("EP_small_10082.off") == 0)
	{
		source_index = 739;
		destination_index = 4444;
	}
	else if (input_dataset.compare("EP_small_20000.off") == 0)
	{
		source_index = 1441;
		destination_index = 8759;
	}
	else if (input_dataset.compare("EP_small_30258.off") == 0)
	{
		source_index = 2141;
		destination_index = 13234;
	}
	else if (input_dataset.compare("EP_small_40328.off") == 0)
	{
		source_index = 2755;
		destination_index = 17693;
	}
	else if (input_dataset.compare("EP_small_50562.off") == 0)
	{
		source_index = 3403;
		destination_index = 22196;
	}
	else
	{
		std::cout << "Dataset not exist!" << std::endl;
		exit(0);
	}

	geodesic::read_mesh_from_file(&input_file[0], points, faces);

	geodesic::Mesh mesh;
	mesh.initialize_mesh_data(points, faces);

	assert(source_index < mesh.vertices().size() && source_index >= 0 &&
		   destination_index < mesh.vertices().size() && destination_index >= 0);

	geodesic::SurfacePoint source = geodesic::SurfacePoint(&mesh.vertices()[source_index]);
	geodesic::SurfacePoint destination = geodesic::SurfacePoint(&mesh.vertices()[destination_index]);

	std::vector<geodesic::SurfacePoint> path;
	std::vector<geodesic::SurfacePoint> path_length_result_path;
	std::vector<geodesic::SurfacePoint> exact_result_path;

	std::vector<geodesic::SurfacePoint> fixed_Steiner_point_result_path;
	std::vector<geodesic::SurfacePoint> fixed_Steiner_point_and_binary_search_snell_law_result_path;
	std::vector<geodesic::SurfacePoint> fixed_Steiner_point_and_effective_weight_snell_law_result_path;

	std::vector<geodesic::SurfacePoint> log_Steiner_point_result_path;
	std::vector<geodesic::SurfacePoint> log_Steiner_point_divide_and_conquer_and_binary_search_snell_law_result_path;
	std::vector<geodesic::SurfacePoint> log_Steiner_point_divide_and_conquer_and_effective_weight_snell_law_result_path;

	int max_loop_num_for_single_endpoint = 30;
	int max_loop_num_for_divide_and_conquer = pow(2, 1);

	assert(Steiner_point_epsilon > 0 && Steiner_point_epsilon <= 1);
	assert(snell_law_epsilon > 0 && snell_law_epsilon <= 1);
	assert(calculate_exact_path == 0 || calculate_exact_path == 1);
	assert(calculate_FSP == 0 || calculate_FSP == 1);

	std::string write_file_header = input_dataset + "\t" +
									std::to_string(mesh.faces().size()) + "\t" +
									std::to_string(std::max(Steiner_point_epsilon, snell_law_epsilon)) + "\t" +
									std::to_string(Steiner_point_epsilon) + "\t" +
									std::to_string(snell_law_epsilon);

	std::ofstream ofs("../output/output.txt", std::ofstream::app);
	ofs << "# dataset\tdatasize\tepsilon\tepsilon_SP\tepsilon_SL\tbuilding_time\tquery_time_step1\tquery_time_step2\tquery_time_total\tmemory_usage_step1\tmemory_usage_step2\tmemory_usage_total\tsnell_law_iteration_count\tdistance_error\tdistance\tedge_sequence_size\n\n";
	ofs.close();

	// simulating the path length
	int estimate_path_length;
	if (!(calculate_exact_path == 0 && calculate_FSP == 0))
	{
		geodesic::GeodesicAlgorithmSubdivision algorithm_path_length(&mesh, 10);
		algorithm_path_length.geodesic(source, destination, path_length_result_path);
		estimate_path_length = path_length_result_path.size();
	}

	// simulating the exact path
	double exact_result_path_distance;
	if (calculate_exact_path == 1)
	{
		std::cout << "\n== Start simulating the exact path ==" << std::endl;
		simulate_exact_path(&mesh, source, destination, path, 0.05, estimate_path_length, 0.05, exact_result_path, exact_result_path_distance);
		std::cout << "== Finished simulating the exact path ==\n"
				  << std::endl;
	}

	if (calculate_FSP == 1)
	{

		std::cout << "== FSP ==" << std::endl;
		fixed_Steiner_point_with_exp_output(&mesh, source, destination, path, write_file_header, Steiner_point_epsilon, estimate_path_length, snell_law_epsilon, exact_result_path_distance, fixed_Steiner_point_result_path);
		std::cout << std::endl;

		std::cout << "== FSP-BSSL ==" << std::endl;
		fixed_Steiner_point_and_binary_search_snell_law(&mesh, source, destination, path, write_file_header, Steiner_point_epsilon, estimate_path_length, snell_law_epsilon, exact_result_path_distance, fixed_Steiner_point_and_binary_search_snell_law_result_path);
		std::cout << std::endl;

		std::cout << "== FSP-EWSL ==" << std::endl;
		fixed_Steiner_point_and_effective_weight_snell_law(&mesh, source, destination, path, write_file_header, Steiner_point_epsilon, estimate_path_length, snell_law_epsilon, exact_result_path_distance, fixed_Steiner_point_and_effective_weight_snell_law_result_path);
		std::cout << std::endl;
	}

	std::cout << "== LSP ==" << std::endl;
	log_Steiner_point_with_exp_output(&mesh, source, destination, path, write_file_header, Steiner_point_epsilon, estimate_path_length, snell_law_epsilon, fixed_Steiner_point_result_path);
	std::cout << std::endl;

	std::cout << "== DLSP-BSSL ==" << std::endl;
	log_Steiner_point_divide_and_conquer_and_binary_search_snell_law(&mesh, source, destination, path, write_file_header, Steiner_point_epsilon, snell_law_epsilon,
																	 max_loop_num_for_divide_and_conquer, max_loop_num_for_single_endpoint, exact_result_path_distance,
																	 log_Steiner_point_divide_and_conquer_and_binary_search_snell_law_result_path);
	std::cout << std::endl;

	std::cout << "== DLSP-EWSL ==" << std::endl;
	log_Steiner_point_divide_and_conquer_and_effective_weight_snell_law(&mesh, source, destination, path, write_file_header, Steiner_point_epsilon, snell_law_epsilon,
																		max_loop_num_for_divide_and_conquer, max_loop_num_for_single_endpoint, exact_result_path_distance,
																		log_Steiner_point_divide_and_conquer_and_effective_weight_snell_law_result_path);
	std::cout << std::endl;

	return 0;
}
